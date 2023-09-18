library(tmvtnorm)
library(dplyr) # For the pipe operator %>%
library(tidyr) # For the function gather()
library(tibble) # For the function rownames_to_column()
library(ggplot2)
library(gridExtra)

truncMVN <- function(n, mu=rep(0, nrow(Sigma)), 
                     Sigma=diag(length(mu)), 
                     a=rep(-Inf, length = length(mu)), 
                     b=rep(Inf, length = length(mu)), 
                     warmup=100, 
                     thinning=1){
  ### n: number of r.v. needed
  ### mu: p-dimensional mean vector
  ### Sigma: p by p Variance-covariance matrix
  ### a: vector of lower bounds; may be -Inf
  ### b: vector of upper bounds; may be Inf
  ### iter: number of iterations used to sample theta
  ### warmup: number of the first most samples discarded
  p = nrow(Sigma)
  Psi = solve(Sigma)
  phi = Psi%*%mu
  d =  max(eigen(Psi)$values) + 1E-5
  A = diag(d, p)-Psi
  L = t(chol(A))
  theta = mu # matrix(mu, nrow=p)
  trace_theta = matrix(NA, nrow=n, ncol=p)
  for (i in (1-warmup):(n*thinning)){
    r <- L%*%rnorm(p) + A%*%theta
    theta <- rtruncnorm(p, a, b, (r+phi)/d, 1/sqrt(d))
    
    ### Burn-in & Thinning
    if(i > 0){
      if (thinning == 1){
        trace_theta[i,] <- theta
      } else if (i %% thinning == 0){
        trace_theta[i %/% thinning,] <- theta
      }
    }
  }
  return(trace_theta)
}



############### 2-d Toy example ###############
par(mfrow=c(1,1))
n_sim = 5
mu = c(0.5, 0.5)
Sigma = matrix(c(1,0.8,0.8,2), 2, 2)
a = c(-1, -Inf)
b = c(0.5, 4)
trace_theta <- matrix(NA, nrow=0, ncol=2)
pb = txtProgressBar(1, n_sim, style=3)
for (i in 1:n_sim){
  trace_theta = rbind(trace_theta, truncMVN(n=1000, mu, Sigma, a, b))
  setTxtProgressBar(pb, i)
}
close(pb)



### Marginal densities
df <- data.frame(theta1=trace_theta[,1], theta2=trace_theta[,2])
df2 <- data.frame(index1=seq(-1,0.5,0.01), d1=dtmvnorm(seq(-1,0.5,0.01), mean=mu, sigma=Sigma, lower=a, upper=b, margin=1))
df3 <- data.frame(index2=seq(-5,4,0.1), d2=dtmvnorm(seq(-5,4,0.1), mean=mu, sigma=Sigma, lower=a, upper=b, margin=2))
gg1 <- ggplot(df, aes(x=theta1)) + geom_density(size=1) + geom_vline(aes(xintercept=-1), color="red", linetype="dashed", size=1) + 
  geom_vline(aes(xintercept=0.5), color="red", linetype="dashed", size=1) +
  geom_line(data=df2, aes(x=index1, y=d1), color="red", size=1) + 
  labs(x=expression(x[1]), y="Density", title=expression(Marginal~density~of~x[1])) +
  theme(plot.title=element_text(size=30, face="bold"), 
        axis.text.x=element_text(size=25), 
        axis.title.x=element_text(size=30), 
        axis.text.y=element_text(size=25), 
        axis.title.y=element_text(size=30))
gg2 <- ggplot(df, aes(x=theta2)) + geom_density() +
  geom_vline(aes(xintercept=4), color="red", linetype="dashed", size=1) + 
  geom_line(data=df3, aes(x=index2, y=d2), color="red", size=1) + 
  labs(x=expression(x[2]), y="Density", title=expression(Marginal~density~of~x[2])) +
  theme(plot.title=element_text(size=30, face="bold"), 
        axis.text.x=element_text(size=25), 
        axis.title.x=element_text(size=30), 
        axis.text.y=element_text(size=25), 
        axis.title.y=element_text(size=30))



### Joint density and contour
contour_lower = c(-5, -5)
contour_upper = c(5, 5)

u <- seq(contour_lower[1], contour_upper[1], 0.05)
v <- seq(contour_lower[2], contour_upper[2], 0.05)

# Generate fake data
den <- function(x, y){
  ### Compute BVN density at (x,y)
  rho = Sigma[1,2] / sqrt(Sigma[1,1]*Sigma[2,2])
  (det(2*pi*Sigma)^(-0.5))*exp(-(((x-mu[1])/Sigma[1,1])^2 + ((y-mu[2])/Sigma[2,2])^2 -2*rho*(x-mu[1])*(y-mu[2])/Sigma[1,1]/Sigma[2,2]) / 2 / (1-rho^2))
}
z = outer(u, v, den)
rownames(z) = u
colnames(z) = v

# Convert data to long format and plot
as.data.frame(z) %>% 
  rownames_to_column(var="row") %>% 
  gather(col, value, -row) %>% 
  mutate(row=as.numeric(row), 
         col=as.numeric(col)) %>% 
  ggplot() + geom_point(data=data.frame(X=c(trace_theta[,1]), Y=c(trace_theta[,2])), aes(x=X, y=Y), size=0.5) + 
  geom_contour(aes(col, row, z=value),bins=20) + 
  labs(x=expression(x[1]), y=expression(x[2]), title="Joint density") +
  theme(plot.title=element_text(size=30), 
        axis.text.x=element_text(size=25), 
        axis.title.x=element_text(size=30), 
        axis.text.y=element_text(size=25), 
        axis.title.y=element_text(size=30)) -> plot_con #+
#theme_classic() -> plot_con
plot_con



grid.arrange(gg1, gg2, plot_con, ncol=3)

