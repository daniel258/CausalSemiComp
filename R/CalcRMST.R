####################################################################################
####################################################################################
# CausalSemiComp
# A function getting a frailty fit and correlation between the frailties
# and return RMSTs
####################################################################################

CalcRMST <- function(rho, tau, n.gamma.vals, n.sample.pers, population, Xnames, data, res, list.out = T, detailed = F)
{
  # Data wrangling to be used for teasing out the baseline curve
  n.X.vals <- nrow(population)
  m.A1 <- mean(data$A==1)
  Xmat <- as.matrix(select(population, Xnames))
  # Calculate exp(X\beta) only once.
  exp.b001 <- exp(Xmat %*% coef(res$fit.list$fit.a0.01))
  exp.b002 <- exp(Xmat %*% coef(res$fit.list$fit.a0.02))
  exp.b012 <- exp(Xmat %*% coef(res$fit.list$fit.a0.12))
  exp.b101 <- exp(Xmat %*% coef(res$fit.list$fit.a1.01))
  exp.b102 <- exp(Xmat %*% coef(res$fit.list$fit.a1.02))
  exp.b112 <- exp(Xmat %*% coef(res$fit.list$fit.a1.12))
  step.A0T1 <- res$H.step.funcs$step.A0T1
  step.A0T2 <- res$H.step.funcs$step.A0T2
  step.A0T12 <- res$H.step.funcs$step.A0T12
  step.A1T1 <- res$H.step.funcs$step.A1T1
  step.A1T2 <- res$H.step.funcs$step.A1T2
  step.A1T12 <- res$H.step.funcs$step.A1T12
  # Replace values larger than tau with tau (later also fix probs)
  A0T1.times <- pmin(knots(step.A0T1), tau) %>% unique
  A0T2.times <- pmin(knots(step.A0T2), tau) %>% unique
  A0T12.times <- pmin(knots(step.A0T12), tau) %>% unique
  A1T1.times <- pmin(knots(step.A1T1), tau) %>% unique
  A1T2.times <- pmin(knots(step.A1T2), tau) %>% unique
  A1T12.times <- pmin(knots(step.A1T12), tau) %>% unique
  # Sample bivariate frailty from gamma distribution
  theta.est <- (1 - m.A1) * res$thetas[1] + m.A1 * res$thetas[2]
  gamma.scale <- theta.est
  gamma.common.shape <- rho/theta.est
  gamma.each.shape <-  (1 - rho)/theta.est
  T1.0.sim <- T1.1.sim <- T2.0.sim <- T2.1.sim <- vector(length = n.X.vals*n.gamma.vals*n.sample.pers)
  # To avoid getting lost: i - which X value j - gamma values wihtin i, k T1,T2 times for specific X, gamma combination
  for(i in 1:n.X.vals)
  {
    #Daniel::CatIndex(i)
    #st <- (i-1) * n.sample.person
    #en <- i * n.sample.person
    gamma.common <- rgamma(n.gamma.vals, shape = gamma.common.shape, scale = gamma.scale)
    gamma0 <- gamma.common + rgamma(n.gamma.vals, shape = gamma.each.shape, scale = gamma.scale)
    gamma1 <- gamma.common + rgamma(n.gamma.vals, shape = gamma.each.shape, scale = gamma.scale)
    for (j in 1:n.gamma.vals)
    {
  #    Daniel::CatIndex(j)
  # For the future: might be better to have linear interpolation between points using approxfun (but f gives global mehtod=constant).
  # For now, only sample from event times
    exp.b001.gamma0.sim <- exp.b001[i] * gamma0[j]
    exp.b002.gamma0.sim <- exp.b002[i] * gamma0[j]
    exp.b012.gamma0.sim <- exp.b012[i] * gamma0[j]
    exp.b101.gamma1.sim <- exp.b101[i] * gamma1[j]
    exp.b102.gamma1.sim <- exp.b102[i] * gamma1[j]
    exp.b112.gamma1.sim <- exp.b112[i] * gamma1[j]
    step.S.a0.01j <- function(t) {exp(-step.A0T1(t) * exp.b001.gamma0.sim)}
    step.S.a0.02j <- function(t) {exp(-step.A0T2(t) * exp.b002.gamma0.sim)}
    step.S.a1.01j <- function(t) {exp(-step.A1T1(t) * exp.b101.gamma1.sim)}
    step.S.a1.02j <- function(t) {exp(-step.A1T2(t) * exp.b102.gamma1.sim)}
    # Assign remaining prob mass to tau
    pr.A0T1 <- -diff(c(1, step.S.a0.01j(A0T1.times)))
    pr.A0T1[length(pr.A0T1)] <- pr.A0T1[length(pr.A0T1)] + 1 - sum(pr.A0T1)
    pr.A0T2 <- -diff(c(1, step.S.a0.02j(A0T2.times)))
    pr.A0T2[length(pr.A0T2)] <- pr.A0T2[length(pr.A0T2)] + 1 - sum(pr.A0T2)
    pr.A1T1 <- -diff(c(1, step.S.a1.01j(A1T1.times)))
    pr.A1T1[length(pr.A1T1)] <- pr.A1T1[length(pr.A1T1)] + 1 - sum(pr.A1T1)
    pr.A1T2 <- -diff(c(1, step.S.a1.02j(A1T2.times)))
    pr.A1T2[length(pr.A1T2)] <- pr.A1T2[length(pr.A1T2)] + 1 - sum(pr.A1T2)
    #### Simulate T1, T2 for n.sample.pers times for each X * gamma values
    T1.0.sim.temp <- sample(A0T1.times, n.sample.pers, replace = T, prob = pr.A0T1)
    T2.0.sim.temp <- sample(A0T2.times, n.sample.pers, replace = T, prob = pr.A0T2)
    T1.1.sim.temp <- sample(A1T1.times, n.sample.pers, replace = T, prob = pr.A1T1)
    T2.1.sim.temp <- sample(A1T2.times, n.sample.pers, replace = T, prob = pr.A1T2)
    ## Redefine T1, and T2:
    # 1. Set T1=Inf for those who were not diseased
    # 2. Resimulate T2 for diseased (or replace with tau)
    disease.0 <- (T1.0.sim.temp < pmin(T2.0.sim.temp, tau))
    disease.1 <- (T1.1.sim.temp < pmin(T2.1.sim.temp, tau))
    T1.0.sim.temp[!disease.0] <- Inf
    T1.1.sim.temp[!disease.1] <- Inf
    T1.0.tau <- T1.0.sim.temp==tau
    T1.1.tau <- T1.1.sim.temp==tau
    T2.0.sim.temp[disease.0 & T1.0.tau] <- tau
    T2.1.sim.temp[disease.1 & T1.1.tau] <- tau
    disease.0.not.tau <- (disease.0 & !T1.0.tau)
    disease.1.not.tau <- (disease.1 & !T1.1.tau)
    s0 <- sum(disease.0.not.tau)
    s1 <- sum(disease.1.not.tau)
    if (s0 > 0)
    {
    for(k in 1:s0)
    {
    #Daniel::CatIndex(k)
    step.S.a0.12k <- function(t) {exp(-(step.A0T12(t) - step.A0T12(T1.0.sim.temp[disease.0.not.tau][k])) * exp.b012.gamma0.sim)}
    A0T12.times.k <- A0T12.times[A0T12.times >= T1.0.sim.temp[disease.0.not.tau][k]]
    pr.A0T12 <- -diff(c(1, step.S.a0.12k(A0T12.times.k)))
    if(length(pr.A0T12)==1) {T2.0.sim.temp[disease.0.not.tau][k] <- tau} else {
    pr.A0T12[length(pr.A0T12)] <- pr.A0T12[length(pr.A0T12)] + 1 - sum(pr.A0T12)
    T2.0.sim.temp[disease.0.not.tau][k] <- sample(A0T12.times.k, 1, replace = T, prob = pr.A0T12)
    }}}
    if (s1 > 0) {
    for(k in 1:s1)
    {
    step.S.a1.12k <- function(t) {exp(-(step.A1T12(t) - step.A1T12(T1.1.sim.temp[disease.1.not.tau][k])) * exp.b112.gamma1.sim)}
    A1T12.times.k <- A1T12.times[A1T12.times >= T1.1.sim.temp[disease.1.not.tau][k]]
    pr.A1T12 <- -diff(c(1, step.S.a1.12k(A1T12.times.k)))
    if(length(pr.A1T12)==1) {T2.1.sim.temp[disease.1.not.tau][k] <- tau} else {
    pr.A1T12[length(pr.A1T12)] <- pr.A1T12[length(pr.A1T12)] + 1 - sum(pr.A1T12)
    T2.1.sim.temp[disease.1.not.tau][k] <- sample(A1T12.times.k, 1, replace = T, prob = pr.A1T12)
    }}}
    st <- (i-1) * n.gamma.vals * n.sample.pers + (j - 1) * n.sample.pers + 1
    en <- (i-1) * n.gamma.vals * n.sample.pers + j * n.sample.pers
    T1.0.sim[st:en] <- T1.0.sim.temp
    T2.0.sim[st:en] <- T2.0.sim.temp
    T1.1.sim[st:en] <- T1.1.sim.temp
    T2.1.sim[st:en] <- T2.1.sim.temp
    }}

  #### Calculate causal effects
  ad <- (T1.0.sim < T2.0.sim) & (T1.1.sim < T2.1.sim)
  nd <- (T1.0.sim >= T2.0.sim) & (T1.1.sim >= T2.1.sim)
  prop.ad <- mean(ad)
  prop.nd <- mean(nd)
  ########################################################
  mean.T2.ad.a1 <- mean(T2.1.sim[ad])
  mean.T2.ad.a0 <- mean(T2.0.sim[ad])
  mean.T2.nd.a1 <- mean(T2.1.sim[nd])
  mean.T2.nd.a0 <- mean(T2.0.sim[nd])
  mean.T1.ad.a1 <- mean(T1.1.sim[ad])
  mean.T1.ad.a0 <- mean(T1.0.sim[ad])
  med.T2.ad.a1 <- median(T2.1.sim[ad])
  med.T2.ad.a0 <- median(T2.0.sim[ad])
  med.T2.nd.a1 <- median(T2.1.sim[nd])
  med.T2.nd.a0 <- median(T2.0.sim[nd])
  med.T1.ad.a1 <- median(T1.1.sim[ad])
  med.T1.ad.a0 <- median(T1.0.sim[ad])
  ATE.T2.ad <- mean.T2.ad.a1 - mean.T2.ad.a0
  ATE.T2.nd <- mean.T2.nd.a1 - mean.T2.nd.a0
  ATE.T1.ad <- mean.T1.ad.a1 - mean.T1.ad.a0
  med.ATE.T2.ad <- med.T2.ad.a1 - med.T2.ad.a0
  med.ATE.T2.nd <- med.T2.nd.a1 - med.T2.nd.a0
  med.ATE.T1.ad <- med.T1.ad.a1 - med.T1.ad.a0
  if (list.out== T) {
    if (detailed == T) {
      ret <- list(mean.T2.ad.a1 = mean.T2.ad.a1, mean.T2.ad.a0 = mean.T2.ad.a0,
                  mean.T1.ad.a1 = mean.T1.ad.a1, mean.T1.ad.a0 = mean.T1.ad.a0,
                  mean.T2.nd.a1 = mean.T2.nd.a1, mean.T2.nd.a0 = mean.T2.nd.a0,
                  med.T2.ad.a1 = med.T1.ad.a1, med.T2.ad.a0 = med.T2.ad.a0,
                  med.T1.ad.a1 = med.T1.ad.a1, med.T1.ad.a0 = med.T1.ad.a0,
                  med.T2.nd.a1 = med.T2.nd.a1, med.T2.nd.a0 = med.T2.nd.a0,
        ATE.T2.ad = ATE.T2.ad, ATE.T2.nd = ATE.T2.nd, ATE.T1.ad = ATE.T1.ad,
                   med.ATE.T2.ad = med.ATE.T2.ad, med.ATE.T2.nd = med.ATE.T2.nd, med.ATE.T1.ad = med.ATE.T1.ad,
        T1.0.sim = T1.0.sim, T1.1.sim = T1.1.sim, T2.0.sim = T2.0.sim, T2.1.sim = T2.1.sim, ad = ad, nd = nd) }
    else {
      ret <- list(ATE.T2.ad = ATE.T2.ad, ATE.T2.nd = ATE.T2.nd, ATE.T1.ad = ATE.T1.ad,
                  med.ATE.T2.ad = med.ATE.T2.ad, med.ATE.T2.nd = med.ATE.T2.nd, med.ATE.T1.ad = med.ATE.T1.ad,
                  prop.ad = prop.ad, prop.nd = prop.nd)
    }
  } else {
    if (detailed == T) {
     ret <- c(mean.T2.ad.a1, mean.T2.ad.a0, mean.T1.ad.a1, mean.T1.ad.a0,
              mean.T2.nd.a1, mean.T2.nd.a0, med.T2.ad.a1, med.T2.ad.a0,
              med.T1.ad.a1, med.T1.ad.a0, med.T2.nd.a1, med.T2.nd.a0,
              ATE.T2.ad, ATE.T2.nd, ATE.T1.ad, med.ATE.T2.ad, med.ATE.T2.nd, med.ATE.T1.ad,
              prop.ad, prop.nd)
     } else {
     ret <- c(ATE.T2.ad, ATE.T2.nd, ATE.T1.ad, med.ATE.T2.ad, med.ATE.T2.nd, med.ATE.T1.ad, prop.ad, prop.nd)
  }}
  return(ret)
}
