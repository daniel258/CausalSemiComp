WrapCausalSCbootSemiParam <- function(data, i = i, tau, rhos, H.times, Xnames, max.iter)
{
  my.data <- data[i, ]
  res <- EMcausalSC(data = my.data, Xnames = Xnames, max.iter = max.iter)
  n.rhos <- length(rhos)
  causal.effects.all.rhos <- vector(length = 6*n.rhos)
  for (j in 1:n.rhos)
  {
    rho <- rhos[j]
    causal.effects <- CalcRMST(rho = rho, tau = rho, n.sample.sim = 100000, data = my.data, Xnames = Xnames, res = res)
    st <- (j-1)*6 +1
    ed <- j*6
    causal.effects.all.rhos[st:ed] <- unlist(causal.effects)
  }
  H.A001 <- res$H.step.funcs$step.A0T1(H.times)
  H.A002 <- res$H.step.funcs$step.A0T2(H.times)
  H.A012 <- res$H.step.funcs$step.A0T12(H.times)
  H.A101 <- res$H.step.funcs$step.A1T1(H.times)
  H.A102 <- res$H.step.funcs$step.A1T2(H.times)
  H.A112 <- res$H.step.funcs$step.A1T12(H.times)
  res.out <- c(res$naive.betas, res$betas, res$thetas, H.A001, H.A002, H.A012, H.A101, H.A102, H.A112,
               causal.effects.all.rhos)
  return(res.out)
}
