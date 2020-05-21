BootEMCausalSC <- function(data, tau, rhos, H.times = NULL, Xnames, max.iter, B)
{
  n.sample <- nrow(data)
  n.rhos <- length(rhos)
  if (!is.null(H.times))
  {
    n.H.times <- length(H.times)
    n.params <- 12 * length(Xnames) + 2 + 6 * n.H.times + 6 * n.rhos # naive betas, betas, thetas, H, causal effects per rho
  } else {
    n.params <- 12 * length(Xnames) + 2 + 6 * n.rhos # naive betas, betas, thetas, causal effects per rho
  }
  all.bs <- matrix(nr = B, nc = n.params)
  for (b in 1:B)
  {
    w <- rexp(n.sample)
    res <- EMcausalSC(data = my.data, Xnames = Xnames, max.iter = max.iter, w = w)
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
    if (!is.null(H.times))
      {
      H.A001 <- res$H.step.funcs$step.A0T1(H.times)
      H.A002 <- res$H.step.funcs$step.A0T2(H.times)
      H.A012 <- res$H.step.funcs$step.A0T12(H.times)
      H.A101 <- res$H.step.funcs$step.A1T1(H.times)
      H.A102 <- res$H.step.funcs$step.A1T2(H.times)
      H.A112 <- res$H.step.funcs$step.A1T12(H.times)
      all.bs[b, ] <- c(res$naive.betas, res$betas, res$thetas, H.A001, H.A002, H.A012, H.A101, H.A102, H.A112,
                 causal.effects.all.rhos)
    } else {
      all.bs[b, ] <- c(res$naive.betas, res$betas, res$thetas, causal.effects.all.rhos)
    }
  }
  all.sd <- apply(all.bs, 2, sd)
  return(all.sd)
}
