WrapEMCausalSCweightedBoot <- function(data, i = i, tau, rhos, H.times, Xnames,
                                       max.iter, init.thetas)
{
  length.out <- 2 + length(Xnames)*6*2 + (length(H.times)*6) # thetas + betas (naive and standard) + H(H.times)
  n.sample <- nrow(data)
  cond <- T
  try <- 1
while (cond==T) { # cond make sure the bootstrap does not fail
  w <- rexp(n.sample)
  res <- tryCatch(EMcausalSC(data = data, Xnames = Xnames, max.iter = max.iter,
                    w = w, init.thetas = init.thetas), error = function(e) {e})
  # if(inherits(res, "error")){
  #   if (try < 3) {try <- try + 1} else {
  #     res.out <- c(rep(NA, length.out), try)
  #     cond <- F
  #   }} else {
  if (any(res$thetas > 20) | any(res$thetas < 0.02))  {
    if (try < 3) {try <- try + 1} else {
      res.out <- c(rep(NA, length.out), try)
      cond <- F
    }} else {
  #n.rhos <- length(rhos)
  #causal.effects.all.rhos <- vector(length = 6*n.rhos)
  # for (j in 1:n.rhos)
  # {
  #   rho <- rhos[j]
  #   causal.effects <- CalcRMST(rho = rho, tau = rho, n.sample.sim = 100000, data = my.data, Xnames = Xnames, res = res)
  #   st <- (j-1)*6 +1
  #   ed <- j*6
  #   causal.effects.all.rhos[st:ed] <- unlist(causal.effects)
  #}
  H.A001 <- res$H.step.funcs$step.A0T1(H.times)
  H.A002 <- res$H.step.funcs$step.A0T2(H.times)
  H.A012 <- res$H.step.funcs$step.A0T12(H.times)
  H.A101 <- res$H.step.funcs$step.A1T1(H.times)
  H.A102 <- res$H.step.funcs$step.A1T2(H.times)
  H.A112 <- res$H.step.funcs$step.A1T12(H.times)
  # res.out <- c(res$naive.betas, res$betas, res$thetas,
  #              H.A001, H.A002, H.A012, H.A101, H.A102, H.A112,
  #              causal.effects.all.rhos)
  res.out <- c(res$naive.betas, res$betas, res$thetas,
               H.A001, H.A002, H.A012, H.A101, H.A102, H.A112, try)
  cond <- F
    }}
  cat(res.out, "\n")
  return(res.out)
}
