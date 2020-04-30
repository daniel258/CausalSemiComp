###########################################################################
# CasualSemiComp

# Data-generating functions for simulations
###########################################################################
# Function that simulate data under frailty Weibull illness-death  model
# with the option to have no disease-protected
# my trick is to simulate more than needed, throw out "dp" and then keep
# only desired sample size
SimDataWeibFrail <- function(n.sample, params, no.protected = T, no.large = T, cens.exp.rate = 0.1,
                             cens.admin = 200)
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
      gamma.out <- cbind(gamma.vec0, gamma.vec1)
      U1.0 <- runif(n.sample.temp)
      U1.1 <- runif(n.sample.temp)
      U2.0 <- runif(n.sample.temp)
      U2.1 <- runif(n.sample.temp)
      T1.0 <- round((-log(U1.0)/(gamma.vec0 * exp(X %*% beta.a0.01) *
                                   base.weib.scale.a0.01^(-base.weib.shape.a0.01)))^
                      (1/base.weib.shape.a0.01), 1)
      T1.1 <- round((-log(U1.1)/(gamma.vec1 * exp(X %*% beta.a1.01) *
                                   base.weib.scale.a1.01^(-base.weib.shape.a1.01)))^
                      (1/base.weib.shape.a1.01), 1)
      T2.0 <- round((-log(U2.0)/(gamma.vec0 * exp(X %*% beta.a0.02) *
                                   base.weib.scale.a0.02^(-base.weib.shape.a0.02)))^
                      (1/base.weib.shape.a0.02), 1)
      T2.1 <- round((-log(U2.1)/(gamma.vec1 * exp(X %*% beta.a1.02) *
                                   base.weib.scale.a1.02^(-base.weib.shape.a1.02)))^
                      (1/base.weib.shape.a1.02), 1)
    } else {
      gamma.vec <- rgamma(n.sample.temp, shape = gamma.shape, scale = gamma.scale)
      #gamma.vec <- rep(1, n.sample.temp)
      gamma.out <- gamma.vec
      U1.0 <- runif(n.sample.temp)
      U1.1 <- runif(n.sample.temp)
      U2.0 <- runif(n.sample.temp)
      U2.1 <- runif(n.sample.temp)
      T1.0 <- round((-log(U1.0)/(gamma.vec * exp(X %*% beta.a0.01) *
                                  base.weib.scale.a0.01^(-base.weib.shape.a0.01)))^
                       (1/base.weib.shape.a0.01), 1)
      T1.1 <- round((-log(U1.1)/(gamma.vec * exp(X %*% beta.a1.01) *
                                  base.weib.scale.a1.01^(-base.weib.shape.a1.01)))^
                      (1/base.weib.shape.a1.01), 1)
      T2.0 <- round((-log(U2.0)/(gamma.vec * exp(X %*% beta.a0.02) *
                                  base.weib.scale.a0.02^(-base.weib.shape.a0.02)))^
                      (1/base.weib.shape.a0.02), 1)
      T2.1 <- round((-log(U2.1)/(gamma.vec * exp(X %*% beta.a1.02) *
                                  base.weib.scale.a1.02^(-base.weib.shape.a1.02)))^
                      (1/base.weib.shape.a1.02), 1)
    }

    #### Simulate Death times for those who were diseased
     for (i in 1:n.sample.temp)
     {
       # if(i==775)
       # {
       #   a <- 3
       #   cat("here!")
       # }
       if (T2.0[i] >= T1.0[i])
       {
         U12.0i <- runif(1)
         T2.0[i] <- round((-log(U12.0i)/(gamma.vec[i] * exp(X[i, ] %*% beta.a0.12) *
                          base.weib.scale.a0.12^(-base.weib.shape.a0.12)) +
                            T1.0[i]^base.weib.shape.a0.12)^(1/base.weib.shape.a0.12), 1)
         if(T2.0[i]==T1.0[i]) {T2.0[i] <- round(T1.0[i] + runif(1, 0.1, 1), 1)}
       }
       if (T2.1[i] >= T1.1[i])
       {
         U12.1i <- runif(1)
         T2.1[i] <- round((-log(U12.1i)/(gamma.vec[i] * exp(X[i, ] %*% beta.a1.12) *
                                           base.weib.scale.a1.12^(-base.weib.shape.a1.12)) +
                             T1.1[i]^base.weib.shape.a1.12)^(1/base.weib.shape.a1.12), 1)
         if(T2.1[i]==T1.1[i]) {T2.1[i] <-  round(T1.1[i] + runif(1,0.1,1), 1)}
       }}
    out.protected <- (T1.0 < T2.0) & (T1.1 > T2.1)
    out.large <- T2.0 > 50 | T2.1 > 50
    if(no.protected==T & no.large==T) { out <- out.protected | out.large}
    if(no.protected==T & no.large==F) { out <- out.protected}
    if(no.protected==F & no.large==T) { out <-  out.large}
    if(any(no.protected, no.large))
    {
      T1.0 <- T1.0[!out]
      T1.1 <- T1.1[!out]
      T2.0 <- T2.0[!out]
      T2.1 <- T2.1[!out]
      X <- X[!out, ]
      gamma.out <- gamma.out[!out]
    }
    n.sample.real <- length(T1.0)
    if(n.sample.real >= n.sample) {cond.sample <- T}
  }
  T1.0 <- T1.0[1:n.sample]
  T1.1 <- T1.1[1:n.sample]
  T2.0 <- T2.0[1:n.sample]
  T2.1 <- T2.1[1:n.sample]
  X <- X[1:n.sample, ]
  gamma.out <- gamma.out[1:n.sample]
  C <- round(rexp(n.sample, rate = cens.exp.rate), 1)
  C[C==T1.0 | C==T2.0 | C==T1.1 | C==T2.1] <- C[C==T1.0 | C==T2.0 | C==T1.1 | C==T2.1] + 0.05
  # Simulate A and obtain observed data
  A <- rbinom(n = n.sample, size = 1, prob = 0.5)
  T1 <- T2 <- delta1 <- delta2 <- vector(length = n.sample)
  T1[A==0] <- pmin(T1.0[A==0], T2.0[A==0], C[A==0], cens.admin)
  T1[A==1] <- pmin(T1.1[A==1], T2.1[A==1], C[A==1], cens.admin)
  T2[A==0] <- pmin(T2.0[A==0], C[A==0], cens.admin)
  T2[A==1] <- pmin(T2.1[A==1], C[A==1], cens.admin)
  delta1[A==0] <-  T1[A==0]==T1.0[A==0]
  delta1[A==1] <-  T1[A==1]==T1.1[A==1]
  delta2[A==0] <-  T2[A==0]==T2.0[A==0]
  delta2[A==1] <-  T2[A==1]==T2.1[A==1]
  T2[delta1==1 & T1==cens.admin] <- cens.admin + 0.05 #avoiding erros in the third move
  list.to.return <- list(T1.0 = T1.0, T1.1 = T1.1, T2.0 = T2.0, T2.1 = T2.1, X = X,
                         T1 = T1, T2 = T2, A = A, C = C, delta1 = delta1, delta2 = delta2,
                         gamma.out = gamma.out)
}
