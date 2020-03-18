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
# without  disease-protected
# my trick is to simulate more than needed, throw out "dp" and then keep
# only desired sample size
SimDataWeibFrailNoDP <- function(n.sample,
                                 base.weib.scale.a0.01, base.weib.scale.a1.01,
                                 base.weib.scale.a0.02, base.weib.scale.a1.02,
                                 base.weib.scale.a0.12, base.weib.scale.a1.12,
                                 base.weib.shape.a0.01, base.weib.shape.a1.01,
                                 base.weib.shape.a0.02, base.weib.shape.a1.02,
                                 base.weib.shape.a0.12, base.weib.shape.a1.12,
                                 theta,
                                 beta.a0.01, beta.a1.01,
                                 beta.a0.02, beta.a1.02,
                                 beta.a0.12, beta.a1.12,
                                 cens.exp.rate = 0.1)
{
  gamma.shape <- 1/theta
  gamma.scale <- theta
  n.sample.temp <- n.sample
  cond.sample <- F
  while (cond.sample==F) {
    n.sample.temp <- n.sample.temp * 4 # is arbitrary
    x1 <- rbinom(n = n.sample.temp, size = 1, prob = 0.5)
    x2 <- rnorm(n.sample.temp)
    X <- cbind(x1, x2)
    gamma.vec <- rgamma(n.sample.temp, shape = gamma.shape, scale = gamma.scale)
    scale.a0.01 <- exp( -log(gamma.vec) - (X %*% beta.a0.01)) * base.weib.scale.a0.01
    scale.a1.01 <- exp( -log(gamma.vec) - (X %*% beta.a1.01)) * base.weib.scale.a1.01
    scale.a0.02 <- exp( -log(gamma.vec) - (X %*% beta.a0.02)) * base.weib.scale.a0.02
    scale.a1.02 <- exp( -log(gamma.vec) - (X %*% beta.a1.02)) * base.weib.scale.a1.02
    scale.a0.12 <- exp( -log(gamma.vec) - (X %*% beta.a0.12)) * base.weib.scale.a0.12
    scale.a1.12 <- exp( -log(gamma.vec) - (X %*% beta.a1.12)) * base.weib.scale.a1.12
    AD.0 <- 65 + round(rweibull(n.sample.temp, shape = base.weib.shape.a0.01, scale = scale.a0.01), 1)
    AD.1 <- 65 + round(rweibull(n.sample.temp, shape = base.weib.shape.a1.01, scale = scale.a1.01), 1)
    Death.0 <- 65 + round(rweibull(n.sample.temp, shape = base.weib.shape.a0.02, scale = scale.a0.02), 1)
    Death.1 <- 65 + round(rweibull(n.sample.temp, shape = base.weib.shape.a1.02, scale = scale.a1.02), 1)
    # while (any(AD.0 > 200)) {
    #   n.sample.temp.AD0 <- sum(AD.0 > 200)
    #   AD.0[AD.0 > 200] <- 65  + round(rweibull(n.sample.temp.AD0,
    #                                            shape = base.weib.shape.a0.02,
    #                                            scale = scale.a0.01[AD.0 > 200]), 1)
    # }
    # while (any(AD.1 > 200)) {
    #   n.sample.temp.AD1 <- sum(AD.1 > 200)
    # AD.1[AD.1 > 200] + round(rweibull(n.sample.temp.AD1,
    #                                   shape = base.weib.shape.a1.02,
    #                                   scale = scale.a1.01[AD.1 > 200]), 1)
    # }
    # while (any(Death.0 > 200)) {
    # n.sample.temp.D0 <- sum(Death.0 > 200)
    #   Death.0[Death.0 > 200] <- 65  + round(rweibull(n.sample.temp.D0,
    #                                                 shape = base.weib.shape.a0.02,
    #                                                 scale = scale.a0.02[Death.0 > 200]), 1)
    # }
    # while (any(Death.1 > 200)) {
    #   n.sample.temp.D1 <- sum(Death.1 > 200)
    #   Death.1[Death.1 > 200] + round(rweibull(n.sample.temp.D1,
    #                                           shape = base.weib.shape.a1.02,
    #                                           scale = scale.a1.02[Death.1 > 200]), 1)
    # }
    #### Simulate Death times for those who were diseased
    for (i in 1:n.sample.temp)
    {
      if (Death.0[i] >= AD.0[i])
       {
         Death.0[i] <- round(CondRweibull(n = 1, shape = base.weib.shape.a0.12,
                                          scale = scale.a0.12[i], t = AD.0[i]), 1)
         # while(Death.0[i] > 200) {
         #       Death.0[i] <- round(CondRweibull(n = 1, shape = base.weib.shape.a0.12,
         #                                        scale = scale.a0.12[i], t = AD.0[i]), 1)
         # }
         if(Death.0[i]==AD.0[i]) {Death.0[i] <- round(AD.0[i] + runif(1,0.1,1), 1)}
       }
       if (Death.1[i] >= AD.1[i])
       {
         Death.1[i] <- round(CondRweibull(n = 1, shape = base.weib.shape.a1.12, scale = scale.a1.12[i], t = AD.1[i]), 1)
         # while(Death.1[i] > 200) {
         #   Death.1[i] <- round(CondRweibull(n = 1, shape = base.weib.shape.a1.12,
         #                                    scale = scale.a1.12[i], t = AD.1[i]), 1)
         # }
         if(Death.1[i]==AD.1[i]) {Death.1[i] <-  round(AD.1[i] + runif(1,0.1,1), 1)}

       }}
    out.protected <- (AD.0 < Death.0) & (AD.1 > Death.1)
    out.large <- Death.0 > 200 | Death.1 > 200
    out <- out.protected | out.large
    AD.0 <- AD.0[!out]
    AD.1 <- AD.1[!out]
    Death.0 <- Death.0[!out]
    Death.1 <- Death.1[!out]
    X <- X[!out, ]
    n.sample.real <- sum(!out)
    if(n.sample.real >= n.sample) {cond.sample <- T}
  }
  AD.0 <- AD.0[1:n.sample]
  AD.1 <- AD.1[1:n.sample]
  Death.0 <- Death.0[1:n.sample]
  Death.1 <- Death.1[1:n.sample]
  X <- X[1:n.sample, ]
  C <- (65 + round(rexp(n.sample, rate = cens.exp.rate), 1))
  # Simulate A and obtain observed data
  A <- rbinom(n = n.sample, size = 1, prob = 0.5)
  T1 <- T2 <- delta1 <- delta2 <- vector(length = n.sample)
  # T1[A==0] <- pmin(AD.0[A==0], Death.0[A==0], C[A==0], 100)
  # T1[A==1] <- pmin(AD.1[A==1], Death.1[A==1], C[A==1], 100)
  # T2[A==0] <- pmin(Death.0[A==0], C[A==0], 100)
  # T2[A==1] <- pmin(Death.1[A==1], C[A==1], 100)
  T1[A==0] <- pmin(AD.0[A==0], Death.0[A==0], C[A==0])
  T1[A==1] <- pmin(AD.1[A==1], Death.1[A==1], C[A==1])
  T2[A==0] <- pmin(Death.0[A==0], C[A==0])
  T2[A==1] <- pmin(Death.1[A==1], C[A==1])
  delta1[A==0] <-  T1[A==0]==AD.0[A==0]
  delta1[A==1] <-  T1[A==1]==AD.1[A==1]
  delta2[A==0] <-  T2[A==0]==Death.0[A==0]
  delta2[A==1] <-  T2[A==1]==Death.1[A==1]
  list.to.return <- list(AD.0 = AD.0, AD.1 = AD.1, Death.0 = Death.0, Death.1 = Death.1,
                         T1 = T1, T2 = T2, A = A, C = C, delta1 = delta1, delta2 = delta2)
}
