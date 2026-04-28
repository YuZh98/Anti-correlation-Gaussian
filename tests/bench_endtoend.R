### End-to-end timing: old vs new l1ball.linreg on the same problem.
###
### Pulls the previous version of SourceCode.R from git history and compares
### running.time and posterior-mean accuracy on a small but realistic problem.

suppressPackageStartupMessages(library(truncnorm))
source("SliceSampler.R")

## --- Materialize the OLD SourceCode.R (one commit before this branch) ---
old_path <- tempfile(fileext = ".R")
git_ref  <- system("git rev-parse HEAD~0 -- SourceCode.R", intern = TRUE)
## We compare against main: the version BEFORE this branch's improvements.
writeLines(system("git show main:SourceCode.R", intern = TRUE), old_path)

## --- Build a problem ---
set.seed(42)
n <- 200; p <- 100
X <- matrix(rnorm(n * p), n, p)
beta_true <- numeric(p); beta_true[1:5] <- c(2, -2, 2, -2, 2)
y <- as.numeric(X %*% beta_true + rnorm(n))
hyp <- list(lam_K0 = 1, aa = 2, bb = 1, cc = 1, dd = 1)

steps  <- 1500
burnin <- 500

## --- OLD ---
old_env <- new.env()
sys.source(old_path, envir = old_env)
environment(old_env$l1ball.linreg) <- old_env
old_env$slice_sampler  <- slice_sampler
old_env$slice_stepout  <- slice_stepout
old_env$slice_shrinkage <- slice_shrinkage
set.seed(7)
sink(tempfile())  # suppress progress bar noise
fit_old <- old_env$l1ball.linreg(X, y, sig20 = 1, tau0 = rep(0.5, p), K0 = 0.5,
                                 hyp = hyp, w = 0.05, m = 10,
                                 steps = steps, burnin = burnin, thin = 1,
                                 init_method = "random")
sink()

## --- NEW ---
source("SourceCode.R")
set.seed(7)
fit_new <- l1ball.linreg(X, y, sig20 = 1, tau0 = rep(0.5, p), K0 = 0.5,
                         hyp = hyp, w = 0.05, slice_m = 10,
                         steps = steps, burnin = burnin, thin = 1,
                         init_method = "random", verbose = FALSE)

cat(sprintf("End-to-end timing (n=%d, p=%d, steps=%d):\n", n, p, steps))
cat(sprintf("  OLD: %.2fs   max|theta_hat - beta_true| = %.3f\n",
            fit_old$running.time, max(abs(fit_old$theta.est - beta_true))))
cat(sprintf("  NEW: %.2fs   max|theta_hat - beta_true| = %.3f\n",
            fit_new$running.time, max(abs(fit_new$theta.est - beta_true))))
cat(sprintf("  speedup: %.2fx\n", fit_old$running.time / fit_new$running.time))
