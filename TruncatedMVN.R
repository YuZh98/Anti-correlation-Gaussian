### Anti-correlation Gaussian sampler for the truncated multivariate normal
### distribution with separable (box) constraints.

if (!requireNamespace("truncnorm", quietly = TRUE)) {
  stop("Package 'truncnorm' is required. Install with install.packages('truncnorm').")
}

truncMVN <- function(n,
                     mu      = rep(0, nrow(Sigma)),
                     Sigma   = diag(length(mu)),
                     a       = rep(-Inf, length(mu)),
                     b       = rep( Inf, length(mu)),
                     warmup  = 100,
                     thinning = 1) {
  ### n:        number of samples to return after warmup/thinning
  ### mu:       p-dimensional mean vector
  ### Sigma:    p x p positive-definite covariance matrix
  ### a, b:     vectors of lower and upper bounds (may be -Inf / Inf)
  ### warmup:   number of warm-up iterations to discard
  ### thinning: keep every thinning-th draw
  ### Returns: an n x p matrix of draws.

  ### ---- Input validation ----
  if (!is.matrix(Sigma)) stop("'Sigma' must be a matrix.")
  p <- nrow(Sigma)
  if (ncol(Sigma) != p)        stop("'Sigma' must be square.")
  if (length(mu) != p)         stop("length(mu) must equal nrow(Sigma).")
  if (length(a)  != p)         stop("length(a) must equal nrow(Sigma).")
  if (length(b)  != p)         stop("length(b) must equal nrow(Sigma).")
  if (any(a >= b))             stop("Need a < b elementwise.")
  if (n      < 1)              stop("'n' must be >= 1.")
  if (warmup < 0)              stop("'warmup' must be >= 0.")
  if (thinning < 1)            stop("'thinning' must be >= 1.")

  Psi <- solve(Sigma)
  phi <- Psi %*% mu
  d   <- max(eigen(Psi, symmetric = TRUE, only.values = TRUE)$values) + 1e-5
  A   <- diag(d, p) - Psi
  L   <- t(chol(A))
  theta <- pmin(pmax(mu, a), b)        # clip the start to the feasible box

  trace_theta <- matrix(NA_real_, nrow = n, ncol = p)
  total_iter  <- warmup + n * thinning
  for (i in seq_len(total_iter)) {
    r     <- L %*% rnorm(p) + A %*% theta
    theta <- truncnorm::rtruncnorm(p, a = a, b = b,
                                   mean = (r + phi) / d, sd = 1 / sqrt(d))
    j <- i - warmup
    if (j > 0) {
      if (thinning == 1) {
        trace_theta[j, ] <- theta
      } else if (j %% thinning == 0) {
        trace_theta[j %/% thinning, ] <- theta
      }
    }
  }
  trace_theta
}


############### 2-d Toy example ###############
if (sys.nframe() == 0L && interactive()) {
  required <- c("tmvtnorm", "dplyr", "tidyr", "tibble", "ggplot2", "gridExtra")
  missing  <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Demo requires packages: ", paste(missing, collapse = ", "))
  }
  library(tmvtnorm)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(gridExtra)

  par(mfrow = c(1, 1))
  n_sim <- 5
  mu    <- c(0.5, 0.5)
  Sigma <- matrix(c(1, 0.8, 0.8, 2), 2, 2)
  a     <- c(-1, -Inf)
  b     <- c(0.5, 4)

  trace_theta <- matrix(NA_real_, nrow = 0, ncol = 2)
  pb <- txtProgressBar(1, n_sim, style = 3)
  for (i in seq_len(n_sim)) {
    trace_theta <- rbind(trace_theta, truncMVN(n = 1000, mu, Sigma, a, b))
    setTxtProgressBar(pb, i)
  }
  close(pb)

  ### Marginal densities (reference is the *truncated* MVN density)
  df  <- data.frame(theta1 = trace_theta[, 1], theta2 = trace_theta[, 2])
  df2 <- data.frame(
    index1 = seq(-1, 0.5, 0.01),
    d1 = dtmvnorm(seq(-1, 0.5, 0.01), mean = mu, sigma = Sigma,
                  lower = a, upper = b, margin = 1)
  )
  df3 <- data.frame(
    index2 = seq(-5, 4, 0.1),
    d2 = dtmvnorm(seq(-5, 4, 0.1), mean = mu, sigma = Sigma,
                  lower = a, upper = b, margin = 2)
  )
  gg1 <- ggplot(df, aes(x = theta1)) + geom_density(size = 1) +
    geom_vline(aes(xintercept = -1), color = "red", linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = 0.5), color = "red", linetype = "dashed", size = 1) +
    geom_line(data = df2, aes(x = index1, y = d1), color = "red", size = 1) +
    labs(x = expression(x[1]), y = "Density",
         title = expression(Marginal ~ density ~ of ~ x[1]))
  gg2 <- ggplot(df, aes(x = theta2)) + geom_density() +
    geom_vline(aes(xintercept = 4), color = "red", linetype = "dashed", size = 1) +
    geom_line(data = df3, aes(x = index2, y = d2), color = "red", size = 1) +
    labs(x = expression(x[2]), y = "Density",
         title = expression(Marginal ~ density ~ of ~ x[2]))

  ### Joint density: use mvtnorm::dmvnorm rather than a hand-rolled formula
  ### (the previous formula divided by variances instead of sds and used the
  ### untruncated density as the reference for truncated samples).
  contour_lower <- c(-5, -5)
  contour_upper <- c( 5,  5)
  u <- seq(contour_lower[1], contour_upper[1], 0.05)
  v <- seq(contour_lower[2], contour_upper[2], 0.05)
  z <- outer(u, v, function(x, y) {
    mvtnorm::dmvnorm(cbind(x, y), mean = mu, sigma = Sigma)
  })
  rownames(z) <- u
  colnames(z) <- v

  plot_con <- as.data.frame(z) %>%
    rownames_to_column(var = "row") %>%
    gather(col, value, -row) %>%
    mutate(row = as.numeric(row), col = as.numeric(col)) %>%
    ggplot() +
    geom_point(data = data.frame(X = trace_theta[, 1], Y = trace_theta[, 2]),
               aes(x = X, y = Y), size = 0.5) +
    geom_contour(aes(col, row, z = value), bins = 20) +
    labs(x = expression(x[1]), y = expression(x[2]), title = "Joint density")

  grid.arrange(gg1, gg2, plot_con, ncol = 3)
}
