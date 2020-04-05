###########################################################################
# CasualSemiComp

# Data-generating functions for simulations
###########################################################################
# A function to simulate Weibull conditional on a value
# Thanks to (http://realizationsinbiostatistics.blogspot.com/2016/05/simulating-weibull-conditional-on-time.html)
CondRweibull <- function(n,shape,scale=1,t=0) {
  if (length(t)!=1 && length(t)!=n) {
    stop("length(t) is not 1 or n")
  }
  return(scale*(-log(runif(n))+(t/scale)^shape)^(1/shape))
}
# Function that simulate data under frailty Weibull illness-death  model
# with the option to have no disease-protected
# my trick is to simulate more than needed, throw out "dp" and then keep
# only desired sample size
SimDataWeibFrail <- function(n.sample, params, no.protected = T, cens.exp.rate = 0.1)
{
  list2env(params, envir = environment())
  # if theta is a scalar same, assume frailty variate is the same. If same dist but indep
  # frail set the same value for bivariate theta vector
  gamma.shape <- 1/theta
  gamma.scale <- theta
  n.sample.temp <- n.sample
  cond.sample <- F
  while (cond.sample==F) {
    n.sample.temp <- n.sample.temp * 4 # is arbitrary
    x1 <- rbinom(n = n.sample.temp, size = 1, prob = 0.5)
    x2 <- rnorm(n.sample.temp)
    X <- cbind(x1, x2)
    if(length(gamma.shape)==2)
    {
      gamma.vec0 <- rgamma(n.sample.temp, shape = gamma.shape[1], scale = gamma.scale[1])
      gamma.vec1 <- rgamma(n.sample.temp, shape = gamma.shape[2], scale = gamma.scale[2])
      scale.a0.01 <- exp( -log(gamma.vec0) - (X %*% beta.a0.01)) * base.weib.scale.a0.01
      scale.a1.01 <- exp( -log(gamma.vec1) - (X %*% beta.a1.01)) * base.weib.scale.a1.01
      scale.a0.02 <- exp( -log(gamma.vec0) - (X %*% beta.a0.02)) * base.weib.scale.a0.02
      scale.a1.02 <- exp( -log(gamma.vec1) - (X %*% beta.a1.02)) * base.weib.scale.a1.02
      scale.a0.12 <- exp( -log(gamma.vec0) - (X %*% beta.a0.12)) * base.weib.scale.a0.12
      scale.a1.12 <- exp( -log(gamma.vec1) - (X %*% beta.a1.12)) * base.weib.scale.a1.12
    } else {
    gamma.vec <- rgamma(n.sample.temp, shape = gamma.shape, scale = gamma.scale)
    scale.a0.01 <- exp( -log(gamma.vec) - (X %*% beta.a0.01)) * base.weib.scale.a0.01
    scale.a1.01 <- exp( -log(gamma.vec) - (X %*% beta.a1.01)) * base.weib.scale.a1.01
    scale.a0.02 <- exp( -log(gamma.vec) - (X %*% beta.a0.02)) * base.weib.scale.a0.02
    scale.a1.02 <- exp( -log(gamma.vec) - (X %*% beta.a1.02)) * base.weib.scale.a1.02
    scale.a0.12 <- exp( -log(gamma.vec) - (X %*% beta.a0.12)) * base.weib.scale.a0.12
    scale.a1.12 <- exp( -log(gamma.vec) - (X %*% beta.a1.12)) * base.weib.scale.a1.12
    }
    T1.0 <- round(rweibull(n.sample.temp, shape = base.weib.shape.a0.01, scale = scale.a0.01), 1)
    T1.1 <- round(rweibull(n.sample.temp, shape = base.weib.shape.a1.01, scale = scale.a1.01), 1)
    T2.0 <- round(rweibull(n.sample.temp, shape = base.weib.shape.a0.02, scale = scale.a0.02), 1)
    T2.1 <- round(rweibull(n.sample.temp, shape = base.weib.shape.a1.02, scale = scale.a1.02), 1)
    #### Simulate Death times for those who were diseased
    for (i in 1:n.sample.temp)
    {
      if (T2.0[i] >= T1.0[i])
       {
         T2.0[i] <- round(CondRweibull(n = 1, shape = base.weib.shape.a0.12,
                                          scale = scale.a0.12[i], t = T1.0[i]), 1)
         if(T2.0[i]==T1.0[i]) {T2.0[i] <- round(T1.0[i] + runif(1,0.1,1), 1)}
       }
       if (T2.1[i] >= T1.1[i])
       {
         T2.1[i] <- round(CondRweibull(n = 1, shape = base.weib.shape.a1.12,
                                       scale = scale.a1.12[i], t = T1.1[i]), 1)
         if(T2.1[i]==T1.1[i]) {T2.1[i] <-  round(T1.1[i] + runif(1,0.1,1), 1)}
       }}
    out.large <- T2.0 > 50 | T2.1 > 50
    if(no.protected==T)
      {
    out.protected <- (T1.0 < T2.0) & (T1.1 > T2.1)
    out <- out.protected | out.large
      }
    else {
      out <- out.large }
    T1.0 <- T1.0[!out]
    T1.1 <- T1.1[!out]
    T2.0 <- T2.0[!out]
    T2.1 <- T2.1[!out]
    X <- X[!out, ]
    n.sample.real <- sum(!out)
    if(n.sample.real >= n.sample) {cond.sample <- T}
  }
  T1.0 <- T1.0[1:n.sample]
  T1.1 <- T1.1[1:n.sample]
  T2.0 <- T2.0[1:n.sample]
  T2.1 <- T2.1[1:n.sample]
  X <- X[1:n.sample, ]
  C <- round(rexp(n.sample, rate = cens.exp.rate), 1)
  # Simulate A and obtain observed data
  A <- rbinom(n = n.sample, size = 1, prob = 0.5)
  T1 <- T2 <- delta1 <- delta2 <- vector(length = n.sample)
  T1[A==0] <- pmin(T1.0[A==0], T2.0[A==0], C[A==0])
  T1[A==1] <- pmin(T1.1[A==1], T2.1[A==1], C[A==1])
  T2[A==0] <- pmin(T2.0[A==0], C[A==0])
  T2[A==1] <- pmin(T2.1[A==1], C[A==1])
  delta1[A==0] <-  T1[A==0]==T1.0[A==0]
  delta1[A==1] <-  T1[A==1]==T1.1[A==1]
  delta2[A==0] <-  T2[A==0]==T2.0[A==0]
  delta2[A==1] <-  T2[A==1]==T2.1[A==1]
  list.to.return <- list(T1.0 = T1.0, T1.1 = T1.1, T2.0 = T2.0, T2.1 = T2.1, X = X,
                         T1 = T1, T2 = T2, A = A, C = C, delta1 = delta1, delta2 = delta2)
}
