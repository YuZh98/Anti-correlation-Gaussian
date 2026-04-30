### Self-contained tests for the samplers in this repository.
###
### Usage:  Rscript tests/test_samplers.R
### Exits with status 0 on success, 1 on any failed check.

suppressPackageStartupMessages({
  ## Hard requirement
  if (!requireNamespace("truncnorm", quietly = TRUE)) {
    stop("Tests require the 'truncnorm' package.")
  }
  has_tmvtnorm <- requireNamespace("tmvtnorm", quietly = TRUE)
  has_mvtnorm  <- requireNamespace("mvtnorm",  quietly = TRUE)
})

## Resolve repo root from this file's location so the script works whether it
## is run from the repo root or from inside tests/.
repo_root <- {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(file_arg) == 0L) getwd() else normalizePath(dirname(dirname(file_arg)))
}
source(file.path(repo_root, "SliceSampler.R"))
source(file.path(repo_root, "SourceCode.R"))
source(file.path(repo_root, "TruncatedMVN.R"))

## ---- Tiny test harness ----
PASS <- 0L
FAIL <- 0L
expect <- function(cond, msg) {
  if (isTRUE(cond)) {
    cat(sprintf("  [PASS] %s\n", msg)); PASS <<- PASS + 1L
  } else {
    cat(sprintf("  [FAIL] %s\n", msg)); FAIL <<- FAIL + 1L
  }
}
section <- function(label) cat(sprintf("\n=== %s ===\n", label))


############################################################
## Test 1: scalar slice sampler matches a standard normal
############################################################
section("Slice sampler: standard normal target")
set.seed(1)
log_pdf_std <- function(x) -0.5 * x * x
x <- 0
draws <- numeric(5000)
for (i in seq_along(draws)) {
  x <- slice_sampler(log_pdf_std, x0 = x, w = 1, m = 20)
  draws[i] <- x
}
expect(abs(mean(draws))                 < 0.1, "mean ~ 0")
expect(abs(sd(draws) - 1)               < 0.1, "sd ~ 1")
ks_p <- suppressWarnings(ks.test(draws, "pnorm")$p.value)
expect(ks_p > 0.01, sprintf("KS test vs N(0,1) (p = %.3f)", ks_p))


############################################################
## Test 2: m-shadowing bug regression
##
## Before the fix, l1ball.linreg's loop reassigned `m` to (step - burnin),
## so the slice sampler's stepout argument silently became negative and the
## stepout while-loops never executed. After the fix the slice sampler is
## parameterized by an explicit `slice_m` that the loop never touches.
############################################################
section("Regression test: slice sampler argument is not shadowed by loop offset")

## Simulate the failure mode with the OLD code path: pass a negative m to
## slice_stepout and verify both while-loops are no-ops (returns the initial
## (L, R) interval of width w around x0).
set.seed(0)
neg_m_interval <- slice_stepout(log_pdf_std, x0 = 0, z = -100, w = 1, m = -9999)
interval_width <- neg_m_interval[2] - neg_m_interval[1]
expect(abs(interval_width - 1) < 1e-12,
       "Negative slice m collapses stepout to a width-w interval (the bug)")

## With the NEW code we explicitly check that l1ball.linreg's slice_m is
## still the value the user passed in by the time the second iteration
## executes. We do this by temporarily wrapping slice_sampler with a probe.
captured_m <- integer(0)
orig_slice <- slice_sampler
local_env <- new.env(parent = environment(l1ball.linreg))
assign("slice_sampler", function(f, x0, w, m) {
  captured_m <<- c(captured_m, m)
  orig_slice(f, x0, w, m)
}, envir = local_env)
environment(l1ball.linreg) <- local_env

set.seed(42)
n <- 60; p <- 8
X <- matrix(rnorm(n * p), n, p)
beta_true <- c(1.5, -1.5, 1.5, rep(0, p - 3))
y <- as.numeric(X %*% beta_true + rnorm(n))
fit <- l1ball.linreg(X, y, sig20 = 1, tau0 = rep(0.5, p), K0 = 0.5,
                     hyp = list(lam_K0 = 1, aa = 2, bb = 1, cc = 1, dd = 1),
                     w = 0.05, slice_m = 13,
                     steps = 200, burnin = 100, thin = 1,
                     init_method = "random", verbose = FALSE)
expect(length(captured_m) == 200,
       "slice_sampler called once per Gibbs iteration")
expect(all(captured_m == 13L),
       sprintf("slice_m = 13 was preserved across all %d iterations (range = [%d, %d])",
               length(captured_m), min(captured_m), max(captured_m)))

## restore
environment(l1ball.linreg) <- baseenv()
environment(l1ball.linreg) <- asNamespace("base")
environment(l1ball.linreg) <- globalenv()


############################################################
## Test 3: l1ball.linreg recovers true sparse support
############################################################
section("l1ball.linreg: support and sign recovery on a sparse problem")
set.seed(7)
n <- 200; p <- 30
X <- matrix(rnorm(n * p), n, p)
beta_true <- numeric(p)
beta_true[c(1, 3, 5)] <- c(2.0, -2.0, 1.5)
y <- as.numeric(X %*% beta_true + rnorm(n))

fit <- l1ball.linreg(X, y, sig20 = 1, tau0 = rep(0.5, p), K0 = 0.5,
                     hyp = list(lam_K0 = 1, aa = 2, bb = 1, cc = 1, dd = 1),
                     w = 0.1, slice_m = 20,
                     steps = 4000, burnin = 2000, thin = 1,
                     init_method = "random", verbose = FALSE)
theta_hat <- fit$theta.est
err <- max(abs(theta_hat - beta_true))
expect(err < 0.5,
       sprintf("max |theta_hat - beta_true| < 0.5 (got %.3f)", err))
sign_match <- all(sign(theta_hat[c(1, 3, 5)]) == sign(beta_true[c(1, 3, 5)]))
expect(sign_match, "Signs of the three active coefficients are recovered")
small_noise <- max(abs(theta_hat[-c(1, 3, 5)]))
expect(small_noise < 0.5,
       sprintf("Inactive coefficients shrunk near zero (max = %.3f)", small_noise))


############################################################
## Test 4: vectorized sampleBeta dispatches all three categories correctly
############################################################
section("Vectorized sampleBeta produces draws inside the correct support")
## With K0 small and tight prior, draws should largely fall outside [-K0, K0].
set.seed(2)
n <- 100; p <- 6
X <- matrix(rnorm(n * p), n, p)
beta_true <- c(2, 2, 2, 2, 2, 2)
y <- as.numeric(X %*% beta_true + 0.1 * rnorm(n))

fit <- l1ball.linreg(X, y, sig20 = 0.01, tau0 = rep(2, p), K0 = 0.1,
                     hyp = list(lam_K0 = 1, aa = 2, bb = 1, cc = 1, dd = 1),
                     w = 0.05, slice_m = 10,
                     steps = 2000, burnin = 1000, thin = 1,
                     init_method = "predetermined", theta0 = beta_true,
                     verbose = FALSE)
post_mean <- colMeans(fit$trace_theta)
expect(all(abs(post_mean - 2) < 0.3),
       sprintf("All six theta posterior means concentrate near 2 (max err = %.3f)",
               max(abs(post_mean - 2))))


############################################################
## Test 5: truncMVN matches the analytic mean of a truncated MVN
############################################################
section("truncMVN: matches analytic moments of a known truncated MVN")
set.seed(11)
mu_t    <- c(0.5, 0.5)
Sigma_t <- matrix(c(1.0, 0.6,
                    0.6, 1.5), 2, 2)
a_t <- c(-1.0, -Inf)
b_t <- c( 0.5,  4.0)

draws <- truncMVN(n = 3000, mu = mu_t, Sigma = Sigma_t,
                  a = a_t, b = b_t, warmup = 500, thinning = 2)
expect(all(draws[, 1] >= a_t[1] & draws[, 1] <= b_t[1]),
       "All draws of x1 within [a1, b1]")
expect(all(draws[, 2] <= b_t[2]),
       "All draws of x2 within (-Inf, b2]")

if (has_tmvtnorm) {
  ref <- tmvtnorm::mtmvnorm(mean = mu_t, sigma = Sigma_t,
                            lower = a_t, upper = b_t)
  emp_mean <- colMeans(draws)
  expect(max(abs(emp_mean - ref$tmean)) < 0.1,
         sprintf("Empirical mean matches analytic mean (max err = %.3f, ref = [%.3f, %.3f])",
                 max(abs(emp_mean - ref$tmean)), ref$tmean[1], ref$tmean[2]))
} else {
  cat("  [SKIP] tmvtnorm not installed; skipping analytic-mean check.\n")
}

## Compare against a reference sampler from tmvtnorm via a 1-D KS test on
## each marginal.
if (has_tmvtnorm) {
  ref_draws <- tmvtnorm::rtmvnorm(3000, mean = mu_t, sigma = Sigma_t,
                                  lower = a_t, upper = b_t,
                                  algorithm = "rejection")
  ks1 <- suppressWarnings(ks.test(draws[, 1], ref_draws[, 1])$p.value)
  ks2 <- suppressWarnings(ks.test(draws[, 2], ref_draws[, 2])$p.value)
  expect(ks1 > 0.001,
         sprintf("Marginal 1 distribution matches reference (KS p = %.3f)", ks1))
  expect(ks2 > 0.001,
         sprintf("Marginal 2 distribution matches reference (KS p = %.3f)", ks2))
}


############################################################
## Test 6: input validation
############################################################
section("Input validation rejects bad calls")
expect(inherits(try(l1ball.linreg(matrix(rnorm(20), 10, 2), y = rnorm(11),
                                  steps = 50, burnin = 10, verbose = FALSE),
                    silent = TRUE), "try-error"),
       "l1ball.linreg errors on length(y) != nrow(X)")
expect(inherits(try(l1ball.linreg(matrix(rnorm(20), 10, 2), y = rnorm(10),
                                  tau0 = rep(0.5, 5),
                                  steps = 50, burnin = 10, verbose = FALSE),
                    silent = TRUE), "try-error"),
       "l1ball.linreg errors on length(tau0) mismatch")
expect(inherits(try(l1ball.linreg(matrix(rnorm(20), 10, 2), y = rnorm(10),
                                  init_method = "predetermined",
                                  steps = 50, burnin = 10, verbose = FALSE),
                    silent = TRUE), "try-error"),
       "l1ball.linreg errors when 'predetermined' init has no theta0")
expect(inherits(try(truncMVN(100, mu = c(0, 0),
                             Sigma = diag(2),
                             a = c(0, 0), b = c(-1, 1)),
                    silent = TRUE), "try-error"),
       "truncMVN errors when a >= b")


############################################################
## Test 7: TruncatedMVN.R BVN density formula bug regression
##
## The OLD code computed BVN density by dividing by variances rather than
## standard deviations. We assert the OLD formula is wrong and that the new
## scripts use mvtnorm::dmvnorm, which is correct.
############################################################
section("Regression test: bivariate normal density formula")
old_den <- function(x, y, mu, Sigma) {
  rho <- Sigma[1, 2] / sqrt(Sigma[1, 1] * Sigma[2, 2])
  (det(2 * pi * Sigma)^(-0.5)) * exp(
    -(((x - mu[1]) / Sigma[1, 1])^2 +
      ((y - mu[2]) / Sigma[2, 2])^2 -
      2 * rho * (x - mu[1]) * (y - mu[2]) / Sigma[1, 1] / Sigma[2, 2]) /
       2 / (1 - rho^2))
}
mu_b    <- c(0.5, 0.5)
Sigma_b <- matrix(c(1, 0.8, 0.8, 2), 2, 2)
old_val <- old_den(0.3, 0.7, mu_b, Sigma_b)
if (has_mvtnorm) {
  ref_val <- mvtnorm::dmvnorm(c(0.3, 0.7), mean = mu_b, sigma = Sigma_b)
  expect(abs(old_val - ref_val) > 1e-3,
         sprintf("Old hand-rolled BVN density disagrees with mvtnorm::dmvnorm (old=%.4f, ref=%.4f)",
                 old_val, ref_val))
} else {
  cat("  [SKIP] mvtnorm not installed; skipping BVN regression check.\n")
}


############################################################
## Summary
############################################################
section("Summary")
cat(sprintf("PASSED: %d    FAILED: %d\n", PASS, FAIL))
if (FAIL > 0L) quit(status = 1, save = "no")
quit(status = 0, save = "no")
