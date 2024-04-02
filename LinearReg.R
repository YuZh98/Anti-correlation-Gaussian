source("~/SliceSampler.R")
source("~/SourceCode.R")


n <- 300
p <- 500
d0 <- 10 # The number of non-zero entries
rho0 <- 0.5 # The correlation of Sigma when x_i\sim N(0,Sigma)
sig2_0 <- 1 # True error variance
c0 <- 3 # The SNR coefficient for |beta_j|

### Initial values
sig20 <- 1
tau0 <- rep(0.5, p)

### Prior parameters
hyp <- list(lam_K0=0.5, aa=5, bb=1, cc=1, dd=1)

### Data Generation
X <- matrix(rnorm(n*p),n)
X0 <- matrix(rnorm(n*p),n)
cord <- c(1:p)
dist <- abs(outer(cord,cord,"-"))
Cor <- (rho0**dist)
X <- X0%*%(chol(Cor))

set.seed(2303)
theta0 <- rep(0, p)
theta0[1:d0] <- c0 * sqrt(sig2_0*log(p)/n) * c(2,-3,2,2,-3,3,-2,3,-2,3) # Ground truth
y0 <- X%*%theta0
y <- sqrt(sig2_0)*rnorm(n) + X%*%theta0


##### Fit the model
fit <- l1ball.linreg(X, y, sig20=sig20, tau0=tau0, hyp=hyp,
                     steps=10000, burnin=2000)
cat("The total running time for the Blocked Gibbs sampler is ", fit$running.time, " seconds.\n")


##### Result analysis
### 1. Estimation
par(mfrow=c(1,1))
if (p<50){
  boxplot(fit$trace_theta, outline=F)
  points(theta0, col="blue", pch=19)
} else{
  boxplot(fit$trace_theta[,1:50], outline=F) # Show first 50 variables
  points(theta0[1:50], col="blue", pch=19)
}

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

### 3. boxplot of ACF
df.ACF = data.frame()
lag.max = 40
for (i in 1:10){
  ACF = acf(fit$trace_theta[,i], lag.max = lag.max, plot=F)
  df.new = data.frame(LAGS=0:(lag.max), acf=ACF$acf[,1,1])
  df.ACF = rbind(df.ACF, df.new)
}
ggplot(df.ACF, aes(x=factor(LAGS), y=acf)) + geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  theme(axis.text.x=element_text(size=25), axis.text.y=element_text(size=25)) +
  scale_x_discrete(breaks = seq(from = 0, to = lag.max, by = 5)) +
  labs(x='', y='')


