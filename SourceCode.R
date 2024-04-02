require(truncnorm)

l1ball.linreg <- function(X, y, sig20=0.1, tau0=rep(0.5,ncol(X)),
                          hyp=list(lam_K0=0.5, aa=5, bb=1, cc=1, dd=1), 
                          steps=3000, burnin=1000, 
                          n=length(y), p=ncol(X)){
  ### Need source('SliceSampler.R') to run this function! ###
  ### X: n by p design matrix
  ### y: n data points
  ### sig20, tau0: initial values for sig2 and tau
  ### hyp: a list of hyperparameters for the prior specification
  ### steps: number of total iterations of the Gibbs sampler
  ### burnin: number of warm-up iterations that are discarded
  ### Output: a list of traces and running time for the Gibbs sampler
  
  ### Getting necessary quantities
  Xy <- t(X)%*%y
  XTX <- crossprod(X)
  M <- t(X)%*%X / sig20
  phi <- Xy / sig20
  lam_p <- max(eigen(XTX)$values)
  E <- diag(lam_p+1E-6, p) - XTX
  LXX <- t(chol(E))
  d = (lam_p+1E-6) / sig20
  
  
  ### Initializing parameters
  sig2 <- sig20
  tau <- tau0
  K0 <- 1
  theta <- rnorm(p)
  beta <- sign(theta)*(abs(theta)+K0)
  
  
  ### Defining auxiliary functions
  softthresholding <- function(beta, threshold){
    (beta-threshold)*(beta>threshold) + (beta+threshold)*(beta<(-threshold))
  }
  
  
  sampleBeta <- function(r){
    temp1 <- d + 1/tau
    temp2 <- sqrt(temp1)
    
    m_0 <- (pnorm(K0/sqrt(tau))-pnorm(-K0/sqrt(tau))) * sqrt(tau)
    m_pos <- pnorm((K0/tau-phi-r)/temp2, lower.tail=F) / temp2
    m_neg <- pnorm((-K0/tau-phi-r)/temp2) / temp2
    
    log_prob_0 <- log(m_0)
    log_prob_pos <- log(m_pos) + (phi+r+d*K0)^2/(2*temp1) - (d*K0^2+2*(phi+r)*K0) / 2
    log_prob_neg <- log(m_neg) + (phi+r-d*K0)^2/(2*temp1) - (d*K0^2-2*(phi+r)*K0) / 2
    
    log_prob_star <- pmax(log_prob_0, log_prob_pos, log_prob_neg)
    sum0 <- exp(log_prob_0-log_prob_star)+exp(log_prob_pos-log_prob_star)+exp(log_prob_neg-log_prob_star)
    sum1 <- exp(log_prob_pos-log_prob_star)+exp(log_prob_neg-log_prob_star)
    
    iszero <- ((log_prob_0-log_prob_star-log(sum0)) > log(runif(p)))
    ispositive <- ((log_prob_pos-log_prob_star-log(sum1)) > log(runif(p)))
    
    beta <- numeric()
    for (j in 1:p){
      if (iszero[j]){
        beta[j] <- rtruncnorm(1, a=-K0, b=K0, mean=0, sd=sqrt(tau[j]))
      } else if (ispositive[j]){
        beta[j] <- rtruncnorm(1, a=K0, b=Inf, mean=(phi[j]+r[j]+d*K0)/temp1[j], sd=1/temp2[j])
      } else {
        beta[j] <- rtruncnorm(1, a=-Inf, b=-K0, mean=(phi[j]+r[j]-d*K0)/temp1[j], sd=1/temp2[j])
      }
    }
    return(beta)
  }
  
  
  sampleR <- function(theta){
    LXX%*%rnorm(p)/sqrt(sig2) + E%*%theta/sig2
  }
  
  
  sampleKappa <- function(kappa0, beta, sig2){
    ### Input "kappa0" and output are both a scalar
    log_f <- function(kappa){
      if (kappa < 0) return(-Inf)
      - sum((y-X%*%softthresholding(beta, kappa))^2) / (2*sig2) - hyp$lam_K0*kappa
    }
    return(slice_sampler(f=log_f, x0=kappa0, w=0.3, m=10))
  }
  
  
  ### Keep track of the trace
  trace_theta <- matrix(0, steps-burnin, p)
  trace_sig2 <- numeric(steps-burnin)
  trace_tau <- matrix(0, steps-burnin, p)
  trace_K0 <- numeric(steps-burnin)
  
  ##### Main
  pb = txtProgressBar(1, steps, style=3)
  running.time <- system.time(for(step in 1:steps){
    ### Gibbs sampler
    r <- sampleR(theta)
    beta <- sampleBeta(r)
    K0 <- sampleKappa(K0, beta, sig2)
    theta <- softthresholding(beta, K0)
    tau <- 1 / rgamma(p, hyp$aa+1/2, rate=hyp$bb+beta^2/2)
    sig2 <- 1 / rgamma(1, hyp$cc+n/2, rate=hyp$dd+sum((y-X%*%theta)^2)/2)
    
    ### Update some quantities
    phi <- Xy / sig2
    d = (lam_p+1E-6)/sig2
    
    ### Burn-in
    m = step-burnin
    if(m>0){
      trace_theta[m,] <- theta
      trace_sig2[m] <- sig2
      trace_tau[m,] <- tau
      trace_K0[m] <- K0
    }
    setTxtProgressBar(pb,step)
  })[3]
  close(pb)
  
  
  return(list(theta.est=colMeans(trace_theta), 
              trace_theta=trace_theta, 
              trace_sig2=trace_sig2,
              trace_tau=trace_tau,
              trace_K0=trace_K0,
              running.time=running.time))
}



