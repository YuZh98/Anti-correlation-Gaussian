### Anti-correlation Gaussian data augmentation Gibbs sampler for the
### L1-ball-prior linear regression model.
###
### Requires SliceSampler.R to be sourced first (provides slice_sampler()).
### Uses truncnorm::rtruncnorm. BoomSpikeSlab is only required when
### init_method = 'SSVS' and is loaded lazily inside that branch.

if (!requireNamespace("truncnorm", quietly = TRUE)) {
  stop("Package 'truncnorm' is required. Install with install.packages('truncnorm').")
}

l1ball.linreg <- function(X, y, sig20 = 0.1, tau0 = rep(0.5, ncol(X)), K0 = 1,
                          hyp = list(lam_K0 = 0.5, aa = 1, bb = 1, cc = 1, dd = 1),
                          w = 0.05, slice_m = 10,
                          steps = 20000, burnin = 10000, thin = 1,
                          init_method = "random", theta0 = NULL,
                          verbose = TRUE) {
  ### X: n by p design matrix
  ### y: n data points
  ### sig20, tau0, K0: initial values for sig2, tau, K0
  ### hyp: a list of hyperparameters for the prior specification
  ### w, slice_m: parameters for the slice sampler (slice_m was named 'm' previously,
  ###             which collided with the loop offset and silently degraded the
  ###             slice sampler's step-out behavior after iteration 1.)
  ### steps: number of total iterations of the Gibbs sampler
  ### burnin: number of warm-up iterations that are discarded
  ### thin: thinning parameter (default is 1, no thinning)
  ### init_method: 'random', 'SSVS', or 'predetermined'
  ### theta0: required if init_method == 'predetermined'
  ### verbose: show a progress bar
  ### Output: a list of traces and running time for the Gibbs sampler

  ### ---- Input validation ----
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.numeric(y)) stop("'y' must be numeric.")
  n <- nrow(X)
  p <- ncol(X)
  if (length(y) != n) stop(sprintf("length(y) = %d but nrow(X) = %d.", length(y), n))
  if (length(tau0) == 1L) tau0 <- rep(tau0, p)
  if (length(tau0) != p) stop(sprintf("length(tau0) must be 1 or %d.", p))
  if (any(tau0 <= 0)) stop("All entries of 'tau0' must be positive.")
  if (sig20 <= 0) stop("'sig20' must be positive.")
  if (K0 < 0) stop("'K0' must be non-negative.")
  if (burnin < 0 || burnin >= steps) stop("Need 0 <= burnin < steps.")
  if (thin < 1 || thin != as.integer(thin)) stop("'thin' must be a positive integer.")
  init_method <- match.arg(init_method, c("random", "SSVS", "predetermined"))

  ### ---- Pre-computation ----
  Xy   <- crossprod(X, y)              # p x 1
  XTX  <- crossprod(X)                 # p x p
  lam_p <- max(eigen(XTX, symmetric = TRUE, only.values = TRUE)$values)
  ridge <- 1e-6
  E    <- diag(lam_p + ridge, p) - XTX
  LXX  <- t(chol(E))

  ### ---- Initialization ----
  sig2 <- sig20
  tau  <- tau0
  phi  <- Xy / sig2
  d    <- (lam_p + ridge) / sig2

  if (init_method == "random") {
    theta <- rnorm(p)
  } else if (init_method == "SSVS") {
    if (!requireNamespace("BoomSpikeSlab", quietly = TRUE)) {
      stop("init_method = 'SSVS' requires the 'BoomSpikeSlab' package.")
    }
    pre_model <- BoomSpikeSlab::lm.spike(y ~ X, niter = 3000)
    theta <- colMeans(pre_model$beta[2000:3000, -1, drop = FALSE])
    sig2  <- mean(pre_model$sigma[2000:3000]^2)
    phi   <- Xy / sig2
    d     <- (lam_p + ridge) / sig2
  } else { # 'predetermined'
    if (is.null(theta0)) stop("theta0 must be provided when init_method = 'predetermined'.")
    if (length(theta0) != p) stop(sprintf("length(theta0) must be %d.", p))
    theta <- theta0
  }
  beta <- sign(theta) * (abs(theta) + K0)

  ### ---- Auxiliary functions ----
  softthresholding <- function(b, threshold) {
    (b - threshold) * (b >  threshold) +
    (b + threshold) * (b < -threshold)
  }

  ### Vectorized sampleBeta. Replaces the original O(p) for-loop with a
  ### single rtruncnorm call per category.
  sampleBeta <- function(r) {
    temp1 <- d + 1 / tau
    temp2 <- sqrt(temp1)
    sqrt_tau <- sqrt(tau)

    ## Component log-weights, computed in the log domain for stability.
    log_m_0   <- log(pnorm(K0 / sqrt_tau) - pnorm(-K0 / sqrt_tau)) + log(sqrt_tau)
    log_m_pos <- pnorm((K0 / tau - phi - r) / temp2, lower.tail = FALSE, log.p = TRUE) -
                 log(temp2)
    log_m_neg <- pnorm((-K0 / tau - phi - r) / temp2, log.p = TRUE) -
                 log(temp2)

    log_prob_0   <- log_m_0
    log_prob_pos <- log_m_pos +
      (phi + r + d * K0)^2 / (2 * temp1) -
      (d * K0^2 + 2 * (phi + r) * K0) / 2
    log_prob_neg <- log_m_neg +
      (phi + r - d * K0)^2 / (2 * temp1) -
      (d * K0^2 - 2 * (phi + r) * K0) / 2

    log_prob_star <- pmax(log_prob_0, log_prob_pos, log_prob_neg)
    sum0 <- exp(log_prob_0   - log_prob_star) +
            exp(log_prob_pos - log_prob_star) +
            exp(log_prob_neg - log_prob_star)
    sum1 <- exp(log_prob_pos - log_prob_star) +
            exp(log_prob_neg - log_prob_star)

    iszero     <- (log_prob_0   - log_prob_star - log(sum0)) > log(runif(p))
    ispositive <- (log_prob_pos - log_prob_star - log(sum1)) > log(runif(p))

    ## Vectorized truncated-normal draws. We pass per-coordinate bounds and
    ## moments and let rtruncnorm dispatch in a single call.
    mean_pos <- (phi + r + d * K0) / temp1
    mean_neg <- (phi + r - d * K0) / temp1
    sd_tail  <- 1 / temp2

    a <- ifelse(iszero, -K0, ifelse(ispositive,  K0, -Inf))
    b <- ifelse(iszero,  K0, ifelse(ispositive,  Inf, -K0))
    mu_vec <- ifelse(iszero, 0, ifelse(ispositive, mean_pos, mean_neg))
    sd_vec <- ifelse(iszero, sqrt_tau, sd_tail)

    truncnorm::rtruncnorm(p, a = a, b = b, mean = mu_vec, sd = sd_vec)
  }

  sampleR <- function(theta) {
    LXX %*% rnorm(p) / sqrt(sig2) + E %*% theta / sig2
  }

  ## sampleKappa is parameterized by w_, m_ explicitly so it does not pick
  ## up shadowed variables from the enclosing loop frame.
  sampleKappa <- function(kappa0, beta, sig2, w_, m_) {
    log_f <- function(kappa) {
      if (kappa < 0) return(-Inf)
      -sum((y - X %*% softthresholding(beta, kappa))^2) / (2 * sig2) -
        hyp$lam_K0 * kappa
    }
    slice_sampler(f = log_f, x0 = kappa0, w = w_, m = m_)
  }

  ### ---- Trace storage ----
  n_samples <- if (thin == 1) (steps - burnin) else ceiling((steps - burnin) / thin)
  trace_theta <- matrix(0, n_samples, p)
  trace_sig2  <- numeric(n_samples)
  trace_tau   <- matrix(0, n_samples, p)
  trace_K0    <- numeric(n_samples)

  ### ---- Main Gibbs loop ----
  if (verbose) pb <- txtProgressBar(1, steps, style = 3)
  running.time <- system.time(for (step in 1:steps) {
    r     <- sampleR(theta)
    beta  <- sampleBeta(r)
    K0    <- sampleKappa(K0, beta, sig2, w_ = w, m_ = slice_m)
    theta <- softthresholding(beta, K0)
    tau   <- 1 / rgamma(p, hyp$aa + 1/2, rate = hyp$bb + beta^2 / 2)

    Xtheta <- X %*% theta                                 # cached
    sig2   <- 1 / rgamma(1, hyp$cc + n/2,
                         rate = hyp$dd + sum((y - Xtheta)^2) / 2)

    phi <- Xy / sig2
    d   <- (lam_p + ridge) / sig2

    ## Burn-in & thinning. Use a non-shadowing local name 'idx'.
    idx <- step - burnin
    if (idx > 0) {
      if (thin == 1) {
        trace_theta[idx, ] <- theta
        trace_sig2[idx]    <- sig2
        trace_tau[idx, ]   <- tau
        trace_K0[idx]      <- K0
      } else if (idx %% thin == 0) {
        k <- idx / thin
        trace_theta[k, ] <- theta
        trace_sig2[k]    <- sig2
        trace_tau[k, ]   <- tau
        trace_K0[k]      <- K0
      }
    }

    if (verbose) setTxtProgressBar(pb, step)
  })[3]
  if (verbose) close(pb)

  list(theta.est    = colMeans(trace_theta),
       trace_theta  = trace_theta,
       trace_sig2   = trace_sig2,
       trace_tau    = trace_tau,
       trace_K0     = trace_K0,
       running.time = running.time)
}
