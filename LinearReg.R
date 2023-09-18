source("~/SliceSampler.R")
source("~/SourceCode.R")

n <- 400
p <- 500
d0 <- 10 # The number of non-zero entries
sig2_0 <- 1E-4 # True error variance


##### Data Generation
X <- matrix(rnorm(n*p),n)
X0 <- matrix(rnorm(n*p),n)
cord <- c(1:p)
dist <- abs(outer(cord,cord,"-"))
Cor <- (0.5**dist)
X <- X0%*%(chol(Cor))
X <- t(t(X)/sqrt(diag(t(X)%*%X)))

set.seed(2303)
theta0 <- rep(0, p)
theta0[1:d0] <- rnorm(d0, 5, 0.5) # Ground truth
y0 <- X%*%theta0
y <- sqrt(sig2_0)*rnorm(n) + X%*%theta0


##### Fit the model
fit <- l1ball.linreg(X, y)
cat("The total running time for the Blocked Gibbs sampler is ", fit$running.time, " seconds.\n")


##### Result analysis
### 1. Estimation
par(mfrow=c(1,1))
plot(theta0, ylab=expression(theta))
points(fit$theta.est, col="red")
legend("topright", legend=c("True theta", "Posterior mean"), col=c("black", "red"), pch=1, cex=1)

### 2. Trace plots (assessing mixing)
N_tol <- length(fit$trace_K0)
## Trace plot of kappa
plot(1:N_tol, fit$trace_K0, "l", main="trace plot of kappa")

## Trace plot of theta
index_to_plot <- c(1:16)
par(mfrow=c(4,4))
for (i in index_to_plot){
  plot(1:N_tol, fit$trace_theta[,i], "l", main=paste("theta", i))
  abline(h=theta0[i], col="red")
}

## Trace plot of sig2
par(mfrow=c(1,1))
plot(1:N_tol, fit$trace_sig2, "l", main="trace plot of sig2")

## Trace plot of tau
par(mfrow=c(4,4))
for (i in index_to_plot){
  plot(1:N_tol, fit$trace_tau[,i], "l", main=paste("tau", i))
}

### 3. ACF plots
for (i in index_to_plot){
  acf(fit$trace_theta[,i])
}


