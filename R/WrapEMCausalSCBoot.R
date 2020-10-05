WrapEMCausalSCBoot <- function(data, i = i, tau = NULL, rhos = NULL, H.times = NULL, Xnames, Lname = NULL,
                               population, max.iter, init.thetas, one.theta = F,
                               n.gamma.vals, n.sample.pers)
{
  library(dplyr)
  library(survival)
  library(CausalSemiComp)
  length.out <- 2 + length(Xnames)*6*2 # thetas + betas (naive and standard)
  if (!is.null(H.times)) {length.out <- length.out + 6*length(H.times)}#+H(H.times)
  if (!is.null(rhos)) {
    n.rhos <- length(rhos)
    length.out <- length.out + 8*n.rhos
    }#+causal effects
  data <- data[i, ]
  n.sample <- nrow(data)
  data$T1 <- data$T1 + runif(n = nrow(data), min = 0, max =  0.0001) # a trick to avoid bootstrap-induced ties
  data$T2 <- data$T2 + runif(n = nrow(data), min = 0, max =  0.0001) # a trick to avoid bootstrap-induced ties
  res <- tryCatch(EMcausalSC(data = data, Xnames = Xnames, Lname = Lname,
                               max.iter = max.iter, init.thetas = init.thetas),
                  error = function(e) {e})

  if(inherits(res, "error"))
    {
      res.out <- rep(NA, length.out)
    } else if (any(res$thetas > 10) | any(res$thetas < 0.02))
      {
        res.out <- rep(NA, length.out)
      }
  else {
      res.out <- c(res$naive.betas, res$betas)
      if(one.theta==T)
      {
        theta.out <-  res$thetas[1] * (1 - res$mean.A) + res$thetas[2]*res$mean.A
      } else {theta.out <- res$thetas}
      res.out <- c(res.out, theta.out)
      if (!is.null(H.times))
      {
      H.A001 <- res$H.step.funcs$step.A0T1(H.times)
      H.A002 <- res$H.step.funcs$step.A0T2(H.times)
      H.A012 <- res$H.step.funcs$step.A0T12(H.times)
      H.A101 <- res$H.step.funcs$step.A1T1(H.times)
      H.A102 <- res$H.step.funcs$step.A1T2(H.times)
      H.A112 <- res$H.step.funcs$step.A1T12(H.times)
      res.out <- c(res.out, H.A001, H.A002, H.A012, H.A101, H.A102, H.A112)
      }
      if (!is.null(rhos)) {
        n.rhos <- length(rhos)
        causal.effects.all.rhos <- vector(length = 8*n.rhos)
        for (j in 1:n.rhos)
        {
          rho <- rhos[j]
          st <- (j-1)*8 +1
          ed <- j*8
          causal.effects.all.rhos[st:ed] <- CalcRMST(rho = rho, tau = tau, n.gamma.vals = n.gamma.vals, n.sample.pers = n.sample.pers,
                                                    population = population, data = data, Xnames = Xnames, res = res,
                                                    list.out = F)
        }
        res.out <- c(res.out, causal.effects.all.rhos)
        }}
  return(res.out)
}
