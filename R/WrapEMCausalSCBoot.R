WrapEMCausalSCBoot <- function(data, i = i, tau = NULL, rhos = NULL, H.times = NULL, Xnames,
                                       max.iter, init.thetas)
{
  library(dplyr)
  library(survival)
  library(CausalSemiComp)
  length.out <- 2 + length(Xnames)*6*2 # thetas + betas (naive and standard)
  if (!is.null(H.times)) {length.out <- length.out + 6*length(H.times)}#+H(H.times)
  if (!is.null(rhos)) {
    n.rhos <- length(rhos)
    length.out <- length.out + 6*n.rhos
    }#+causal effects
  data <- data[i, ]
  res <- tryCatch(EMcausalSC(data = data, Xnames = Xnames, max.iter = max.iter,
                             init.thetas = init.thetas), error = function(e) {e})
  if(inherits(res, "error"))
    {
      res.out <- rep(NA, length.out)
    } else if (any(res$thetas > 10) | any(res$thetas < 0.02))
      {
        res.out <- rep(NA, length.out)
      }
  else {
      res.out <- c(res$naive.betas, res$betas, res$thetas)
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
        causal.effects.all.rhos <- vector(length = 6*n.rhos)
        for (j in 1:n.rhos)
        {
          rho <- rhos[j]
          causal.effects <- CalcRMST(rho = rho, tau = rho, n.sample.sim = 100000, data = my.data, Xnames = Xnames, res = res)
          st <- (j-1)*6 +1
          ed <- j*6
          causal.effects.all.rhos[st:ed] <- unlist(causal.effects)
          res.out <- c(res.out, causal.effects.all.rhos)
        }
        }}
  return(res.out)
}
