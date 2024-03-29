#' @title Functions simulating data under frailty Weibull illness-death model
#  with the option to have no disease-protected
#' @description Acceptance-rejection sampling is used in case the order preservation assumption is imposed.
#' @n.sample Desired sample size
#' @params A list of Weibull illness-death model; see \code{Weibull illness-death model}
#' @author Daniel Nevo
#' @importFrom stats rgamma rbinom runif rexp rnorm
#' @export
SimDataWeibFrail <- function(n.sample, params, no.protected = T, no.large = T, cens.exp.rate = 0.1,
                             round.times = T, cens.admin = 100, X = NULL)
{
  if (!is.null(X) & (no.large | no.protected)) stop("Can't specify X and ask no large or no protected")
  list2env(params, envir = environment())
  gamma.scale <- params$theta
  gamma.common.shape <- params$rho/params$theta
  gamma.each.shape <-  (1 - params$rho)/params$theta
  n.sample.temp <- n.sample
  cond.sample <- F
  while (cond.sample==F) {
    if(no.protected | no.large) {
    n.sample.temp <- n.sample.temp * 4 # is arbitrary - this is a trick for the case no.protected or no.larger is TRUE
    }
    if (is.null(X))
    {
      x1 <- rbinom(n = n.sample.temp, size = 1, prob = 0.5)
      x2 <- rnorm(n.sample.temp)
      X <- cbind(x1, x2)
    } else {
      n.times <- ceiling(n.sample.temp/nrow(X))
      X <- do.call(rbind, replicate(n.times, X, simplify=F)) %>% as.matrix
      n.sample.temp <- nrow(X)
    }
    gamma.common <- rgamma(n.sample.temp, shape = gamma.common.shape, scale = gamma.scale)
    gamma0 <- gamma.common + rgamma(n.sample.temp, shape = gamma.each.shape, scale = gamma.scale)
    gamma1 <- gamma.common + rgamma(n.sample.temp, shape = gamma.each.shape, scale = gamma.scale)
    gamma.out <- cbind(gamma0, gamma1)
    U1.0 <- runif(n.sample.temp)
    U1.1 <- runif(n.sample.temp)
    U2.0 <- runif(n.sample.temp)
    U2.1 <- runif(n.sample.temp)
    U12.0 <- runif(n.sample.temp)
    U12.1 <- runif(n.sample.temp)
    expb001.gamma <- gamma0 * exp(X %*% params$beta.a0.01)
    expb002.gamma <- gamma0 * exp(X %*% params$beta.a0.02)
    expb012.gamma <- gamma0 * exp(X %*% params$beta.a0.12)
    expb101.gamma <- gamma1 * exp(X %*% params$beta.a1.01)
    expb102.gamma <- gamma1 * exp(X %*% params$beta.a1.02)
    expb112.gamma <- gamma1 * exp(X %*% params$beta.a1.12)
    T1.0 <- (-log(U1.0)/(expb001.gamma * params$base.weib.scale.a0.01^(-params$base.weib.shape.a0.01)))^
                              (1/params$base.weib.shape.a0.01)
    T1.1 <- (-log(U1.1)/(expb101.gamma * params$base.weib.scale.a1.01^(-params$base.weib.shape.a1.01)))^
                      (1/params$base.weib.shape.a1.01)
    T2.0 <- (-log(U2.0)/(expb002.gamma * params$base.weib.scale.a0.02^(-params$base.weib.shape.a0.02)))^
                      (1/params$base.weib.shape.a0.02)
    T2.1 <- (-log(U2.1)/(expb102.gamma * params$base.weib.scale.a1.02^(-params$base.weib.shape.a1.02)))^
                      (1/params$base.weib.shape.a1.02)
    if (round.times==T) {
      T1.0 <- round(T1.0, 1)
      T1.1 <- round(T1.1, 1)
      T2.0 <- round(T2.0, 1)
      T2.1 <- round(T2.1, 1)
    }
    #### Re-simulate death times for those who were diseased
    for (i in 1:n.sample.temp)
    {
      if (T2.0[i] >= T1.0[i])
      {
        U12.0i <- U12.0[i]
        T2.0[i] <- (-log(U12.0i)/(expb012.gamma[i] * params$base.weib.scale.a0.12^(-params$base.weib.shape.a0.12)) +
                     T1.0[i]^params$base.weib.shape.a0.12)^(1/params$base.weib.shape.a0.12)
        if (round.times==T) {
          T2.0[i] <- round(T2.0[i], 1)
        if(T2.0[i]==T1.0[i]) {T2.0[i] <- round(T1.0[i] + runif(1, 0.1, 1), 1)}
        }
       }
      if (T2.1[i] >= T1.1[i])
      {
        U12.1i <- U12.1[i]
        T2.1[i] <- (-log(U12.1i)/(expb112.gamma[i] * params$base.weib.scale.a1.12^(-params$base.weib.shape.a1.12)) +
                             T1.1[i]^params$base.weib.shape.a1.12)^(1/params$base.weib.shape.a1.12)
        if (round.times==T) {
          T2.1[i] <- round(T2.1[i], 1)
          if(T2.1[i]==T1.1[i]) {T2.1[i] <- round(T1.1[i] + runif(1, 0.1, 1), 1)}
        }
       }}
    out.protected <- (T1.0 < T2.0) & (T1.1 > T2.1)
    out.large <- T2.0 > 50 | T2.1 > 50
    if(no.protected==T & no.large==T) {out <- out.protected | out.large}
    if(no.protected==T & no.large==F) {out <- out.protected}
    if(no.protected==F & no.large==T) {out <-  out.large}
    if(any(no.protected, no.large))
    {
      T1.0 <- T1.0[!out]
      T1.1 <- T1.1[!out]
      T2.0 <- T2.0[!out]
      T2.1 <- T2.1[!out]
      X <- X[!out, ]
      gamma.out <- gamma.out[!out, ]
    }
    n.sample.real <- length(T1.0)
    if(n.sample.real >= n.sample) {cond.sample <- T}
  }
  T1.0 <- T1.0[1:n.sample]
  T1.1 <- T1.1[1:n.sample]
  T2.0 <- T2.0[1:n.sample]
  T2.1 <- T2.1[1:n.sample]
  # Fix zeros
  if (round.times==T) {
  T2.0[T1.0==0] <- pmax(T2.0[T1.0==0], 0.1)
  T2.1[T1.1==0] <- pmax(T2.1[T1.1==0], 0.1)
  T1.0[T1.0==0] <- 0.1
  T1.1[T1.1==0] <- 0.1
  T2.0[T2.0==0] <- 0.1
  T2.1[T2.1==0] <- 0.1
  }
  X <- X[1:n.sample, ]
  gamma.out <- gamma.out[1:n.sample, ]
  C <- rexp(n.sample, rate = cens.exp.rate)
  if (round.times==T) {
  C <- round(C, 1)
  C[C==T1.0 | C==T2.0 | C==T1.1 | C==T2.1] <- C[C==T1.0 | C==T2.0 | C==T1.1 | C==T2.1] + 0.05
  }
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


# The function below was used for internal checkings and should not be used without careful inspection.
SimDataWeibFrailL <- function(n.sample, params, no.protected = T, no.large = T, cens.exp.rate = 0.1,
                             round.times = T, cens.admin = 100, X = NULL, Lmax = NULL)
{
  if (!is.null(X) & (no.large | no.protected)) stop("Can't specify X and ask no large or no protected")
  list2env(params, envir = environment())
  gamma.scale <- params$theta
  gamma.common.shape <- params$rho/params$theta
  gamma.each.shape <-  (1 - params$rho)/params$theta
  n.sample.temp <- n.sample
  cond.sample <- F
  while (cond.sample==F) {
    if(no.protected | no.large) {
      n.sample.temp <- n.sample.temp * 4 # is arbitrary - this is a trick for the case no.protected or no.larger is TRUE
    }
    if (is.null(X))
    {
      x1 <- rbinom(n = n.sample.temp, size = 1, prob = 0.5)
      x2 <- rnorm(n.sample.temp)
      X <- cbind(x1, x2)
    } else {
      n.times <- ceiling(n.sample.temp/nrow(X))
      X <- do.call(rbind, replicate(n.times, X, simplify=F)) %>% as.matrix
      n.sample.temp <- nrow(X)
    }
    gamma.common <- rgamma(n.sample.temp, shape = gamma.common.shape, scale = gamma.scale)
    gamma0 <- gamma.common + rgamma(n.sample.temp, shape = gamma.each.shape, scale = gamma.scale)
    gamma1 <- gamma.common + rgamma(n.sample.temp, shape = gamma.each.shape, scale = gamma.scale)
    gamma.out <- cbind(gamma0, gamma1)
    U1.0 <- runif(n.sample.temp)
    U1.1 <- runif(n.sample.temp)
    U2.0 <- runif(n.sample.temp)
    U2.1 <- runif(n.sample.temp)
    U12.0 <- runif(n.sample.temp)
    U12.1 <- runif(n.sample.temp)
    expb001.gamma <- gamma0 * exp(X %*% params$beta.a0.01)
    expb002.gamma <- gamma0 * exp(X %*% params$beta.a0.02)
    expb012.gamma <- gamma0 * exp(X %*% params$beta.a0.12)
    expb101.gamma <- gamma1 * exp(X %*% params$beta.a1.01)
    expb102.gamma <- gamma1 * exp(X %*% params$beta.a1.02)
    expb112.gamma <- gamma1 * exp(X %*% params$beta.a1.12)
    T1.0 <- (-log(U1.0)/(expb001.gamma * params$base.weib.scale.a0.01^(-params$base.weib.shape.a0.01)))^
      (1/params$base.weib.shape.a0.01)
    T1.1 <- (-log(U1.1)/(expb101.gamma * params$base.weib.scale.a1.01^(-params$base.weib.shape.a1.01)))^
      (1/params$base.weib.shape.a1.01)
    T2.0 <- (-log(U2.0)/(expb002.gamma * params$base.weib.scale.a0.02^(-params$base.weib.shape.a0.02)))^
      (1/params$base.weib.shape.a0.02)
    T2.1 <- (-log(U2.1)/(expb102.gamma * params$base.weib.scale.a1.02^(-params$base.weib.shape.a1.02)))^
      (1/params$base.weib.shape.a1.02)
    if (round.times==T) {
      T1.0 <- round(T1.0, 1)
      T1.1 <- round(T1.1, 1)
      T2.0 <- round(T2.0, 1)
      T2.1 <- round(T2.1, 1)
    }
    #### Re-simulate death times for those who were diseased
    for (i in 1:n.sample.temp)
    {
      if (T2.0[i] >= T1.0[i])
      {
        U12.0i <- U12.0[i]
        T2.0[i] <- (-log(U12.0i)/(expb012.gamma[i] * params$base.weib.scale.a0.12^(-params$base.weib.shape.a0.12)) +
                      T1.0[i]^params$base.weib.shape.a0.12)^(1/params$base.weib.shape.a0.12)
        if (round.times==T) {
          T2.0[i] <- round(T2.0[i], 1)
          if(T2.0[i]==T1.0[i]) {T2.0[i] <- round(T1.0[i] + runif(1, 0.1, 1), 1)}
        }
      }
      if (T2.1[i] >= T1.1[i])
      {
        U12.1i <- U12.1[i]
        T2.1[i] <- (-log(U12.1i)/(expb112.gamma[i] * params$base.weib.scale.a1.12^(-params$base.weib.shape.a1.12)) +
                      T1.1[i]^params$base.weib.shape.a1.12)^(1/params$base.weib.shape.a1.12)
        if (round.times==T) {
          T2.1[i] <- round(T2.1[i], 1)
          if(T2.1[i]==T1.1[i]) {T2.1[i] <- round(T1.1[i] + runif(1, 0.1, 1), 1)}
        }
      }}
    out.protected <- (T1.0 < T2.0) & (T1.1 > T2.1)
    out.large <- T2.0 > 50 | T2.1 > 50
    if(no.protected==T & no.large==T) {out <- out.protected | out.large}
    if(no.protected==T & no.large==F) {out <- out.protected}
    if(no.protected==F & no.large==T) {out <-  out.large}
    if(any(no.protected, no.large))
    {
      T1.0 <- T1.0[!out]
      T1.1 <- T1.1[!out]
      T2.0 <- T2.0[!out]
      T2.1 <- T2.1[!out]
      X <- X[!out, ]
      gamma.out <- gamma.out[!out, ]
    }
    n.sample.real <- length(T1.0)
    if(n.sample.real >= n.sample) {cond.sample <- T}
  }
  T1.0 <- T1.0[1:n.sample]
  T1.1 <- T1.1[1:n.sample]
  T2.0 <- T2.0[1:n.sample]
  T2.1 <- T2.1[1:n.sample]
  # Fix zeros
  if (round.times==T) {
    T2.0[T1.0==0] <- pmax(T2.0[T1.0==0], 0.1)
    T2.1[T1.1==0] <- pmax(T2.1[T1.1==0], 0.1)
    T1.0[T1.0==0] <- 0.1
    T1.1[T1.1==0] <- 0.1
    T2.0[T2.0==0] <- 0.1
    T2.1[T2.1==0] <- 0.1
  }
  X <- X[1:n.sample, ]
  gamma.out <- gamma.out[1:n.sample, ]
  C <- rexp(n.sample, rate = cens.exp.rate)
  if (round.times==T) {
    C <- round(C, 1)
    C[C==T1.0 | C==T2.0 | C==T1.1 | C==T2.1] <- C[C==T1.0 | C==T2.0 | C==T1.1 | C==T2.1] + 0.05
  }
  # Simulate L, A and obtain observed data
  if (!is.na(Lmax))
  {
    L <- runif(n.sample, 0, Lmax)
  }
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
  not.truncated <- (L < T1) & (L < C)
  L <- L[not.truncated]
  A <- A[not.truncated]
  X <- X[not.truncated, ]
  T1 <- T1[not.truncated]
  T2 <- T2[not.truncated]
  delta1 <- delta1[not.truncated]
  delta2 <- delta2[not.truncated]
  list.to.return <- list(T1.0 = T1.0, T1.1 = T1.1, T2.0 = T2.0, T2.1 = T2.1, X = X, L = L,
                         T1 = T1, T2 = T2, A = A, C = C, delta1 = delta1, delta2 = delta2,
                         gamma.out = gamma.out)
}
