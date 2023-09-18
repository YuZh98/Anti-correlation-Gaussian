##### Scalar slice sampler
slice_stepout <- function(f, x0, z, w=0.5, m=10){
  ### f = logarithm of the function proportional to the density
  ### x0 = the current point
  ### z = the vertical level defining the slice
  ### w = estimate of the typical size of a slice
  ### m = integer limiting the size of a slice to mw
  ### Output: (L,R) = the interval found
  L <- x0 - w*runif(1)
  R <- L + w
  J <- floor(m*runif(1))
  K <- m - 1 - J
  while (J>0 & z<f(L)){
    L <- L - w
    J <- J - 1
  }
  while (K>0 & z<f(R)){
    R <- R + w
    K <- K - 1
  }
  return(c(L, R))
}


slice_shrinkage <- function(f, x0, z, I){
  ### f = logarithm of the function proportional to the density
  ### x0 = the current point
  ### z = the vertical level defining the slice
  ### I = (L, R) = the interval to sample from
  ### Output: x1 = the new point
  Lbar <- I[1]
  Rbar <- I[2]
  repeat{
    x1 <- Lbar + runif(1)*(Rbar-Lbar)
    if (z<f(x1)){
      return(x1)
    }
    if (x1<x0){
      Lbar <- x1
    } else{
      Rbar <- x1
    }
  }
}


slice_sampler <- function(f, x0, w=0.5, m=10){
  ### f = logarithm of the function proportional to the density
  ### x0 = the current point
  ### w = estimate of the typical size of a slice
  ### m = integer limiting the size of a slice to mw
  ### Output: x1 = the new point
  z <- f(x0) - rexp(1)
  I <- slice_stepout(f, x0, z, w, m)
  return(slice_shrinkage(f, x0, z, I))
}






##### Multi-dimensional slice sampler
dir_slice_stepout <- function(f, x0, dir, z, w=0.5, m=10){
  ### f = logarithm of the function proportional to the density
  ### x0 = the current point
  ### dir = the updating direction
  ### z = the vertical level defining the slice
  ### w = estimate of the typical size of a slice
  ### m = integer limiting the size of a slice to mw
  ### Output: (L,R) = the intervals found
  L <- x0 - w*runif(1)*dir
  R <- L + w*dir
  J <- floor(m*runif(1))
  K <- m - 1 - J
  while (J>0 & z<f(L)){
    L <- L - w*dir
    J <- J - 1
  }
  while (K>0 & z<f(R)){
    R <- R + w*dir
    K <- K - 1
  }
  return(cbind(L, R))
}


dir_slice_shrinkage <- function(f, x0, dir, z, I){
  ### f = logarithm of the function proportional to the density
  ### x0 = the current point
  ### dir = the updating direction
  ### z = the vertical level defining the slice
  ### I = (L, R) = the intervals to sample from
  ### Output: x1 = the new point
  Lbar <- c(I[,1])
  Rbar <- c(I[,2])
  repeat{
    x1 <- Lbar + runif(1)*(Rbar-Lbar)
    if (z<f(x1)){
      return(x1)
    }
    if (sum((x1-x0)*dir)<0){
      Lbar <- x1
    } else{
      Rbar <- x1
    }
  }
}


dir_slice_sampler <- function(f, x0, dir, w=0.5, m=10){
  ### f = logarithm of the function proportional to the density
  ### x0 = the current point
  ### dir = the updating direction
  ### w = estimate of the typical size of a slice
  ### m = integer limiting the size of a slice to mw
  ### Output: x1 = the new point
  z <- f(x0) - rexp(1)
  I <- dir_slice_stepout(f, x0, dir, z, w, m)
  return(dir_slice_shrinkage(f, x0, dir, z, I))
}









