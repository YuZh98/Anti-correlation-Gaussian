### Micro-benchmark: vectorized vs scalar-loop sampleBeta.
###
### Reproduces the inner loop structure of l1ball.linreg::sampleBeta and
### times both implementations on the same inputs. Useful as a sanity check
### that the vectorized version is genuinely faster.

suppressPackageStartupMessages(library(truncnorm))
set.seed(0)

p     <- 500
tau   <- runif(p, 0.5, 2)
phi   <- rnorm(p)
r     <- rnorm(p)
d     <- 3
K0    <- 0.5

temp1 <- d + 1 / tau
temp2 <- sqrt(temp1)
sqrt_tau <- sqrt(tau)

m_0   <- (pnorm(K0 / sqrt_tau) - pnorm(-K0 / sqrt_tau)) * sqrt_tau
m_pos <- pnorm((K0 / tau - phi - r) / temp2, lower.tail = FALSE) / temp2
m_neg <- pnorm((-K0 / tau - phi - r) / temp2) / temp2
log_prob_0   <- log(m_0)
log_prob_pos <- log(m_pos) +
  (phi + r + d * K0)^2 / (2 * temp1) -
  (d * K0^2 + 2 * (phi + r) * K0) / 2
log_prob_neg <- log(m_neg) +
  (phi + r - d * K0)^2 / (2 * temp1) -
  (d * K0^2 - 2 * (phi + r) * K0) / 2
log_prob_star <- pmax(log_prob_0, log_prob_pos, log_prob_neg)
sum0 <- exp(log_prob_0   - log_prob_star) +
        exp(log_prob_pos - log_prob_star) +
        exp(log_prob_neg - log_prob_star)
sum1 <- exp(log_prob_pos - log_prob_star) +
        exp(log_prob_neg - log_prob_star)

set.seed(1)
iszero     <- (log_prob_0   - log_prob_star - log(sum0)) > log(runif(p))
ispositive <- (log_prob_pos - log_prob_star - log(sum1)) > log(runif(p))

scalar_loop <- function() {
  beta <- numeric(p)
  for (j in 1:p) {
    if (iszero[j]) {
      beta[j] <- rtruncnorm(1, a = -K0, b = K0, mean = 0, sd = sqrt(tau[j]))
    } else if (ispositive[j]) {
      beta[j] <- rtruncnorm(1, a = K0, b = Inf,
                            mean = (phi[j] + r[j] + d * K0) / temp1[j],
                            sd   = 1 / temp2[j])
    } else {
      beta[j] <- rtruncnorm(1, a = -Inf, b = -K0,
                            mean = (phi[j] + r[j] - d * K0) / temp1[j],
                            sd   = 1 / temp2[j])
    }
  }
  beta
}

vectorized <- function() {
  mean_pos <- (phi + r + d * K0) / temp1
  mean_neg <- (phi + r - d * K0) / temp1
  sd_tail  <- 1 / temp2
  a <- ifelse(iszero, -K0, ifelse(ispositive,  K0, -Inf))
  b <- ifelse(iszero,  K0, ifelse(ispositive,  Inf, -K0))
  mu_vec <- ifelse(iszero, 0, ifelse(ispositive, mean_pos, mean_neg))
  sd_vec <- ifelse(iszero, sqrt_tau, sd_tail)
  rtruncnorm(p, a = a, b = b, mean = mu_vec, sd = sd_vec)
}

reps <- 200
t_scalar <- system.time(for (i in seq_len(reps)) { set.seed(123); scalar_loop() })
t_vec    <- system.time(for (i in seq_len(reps)) { set.seed(123); vectorized()  })

cat(sprintf("p = %d, reps = %d\n", p, reps))
cat(sprintf("scalar loop:  elapsed = %.3fs\n",  t_scalar["elapsed"]))
cat(sprintf("vectorized:   elapsed = %.3fs\n",  t_vec["elapsed"]))
cat(sprintf("speedup:      %.2fx\n", t_scalar["elapsed"] / t_vec["elapsed"]))
