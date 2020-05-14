####################################################################################
####################################################################################
# CausalSemiComp
# A function getting a frailty fit and correlation between the frailties
# and return RMSTs
####################################################################################

CalcRMST <- function(rho, tau, n.sample.sim, data, Xnames, res)
{

  if(n.sample.sim < 10001) {warning(paste0("Hard to believe that n.sample = ", n.sample.sim,
  "is sufficient to approximate infinite population"))}
  # Data wrangling to be used for teasing out the baseline curve
  data.predictA0 <- data %>% dplyr::filter(A==0)
  data.predictA1 <- data %>% dplyr::filter(A==1)
  data.predictA0T12 <- data %>% dplyr::filter(A==0 & delta1==1)
  data.predictA1T12 <- data %>% dplyr::filter(A==1 & delta1==1)
  XmatA0 <- as.matrix(select(data.predictA0, Xnames))
  XmatA1 <- as.matrix(select(data.predictA1, Xnames))
  XmatA0T12 <- as.matrix(select(data.predictA0T12, Xnames))
  XmatA1T12 <- as.matrix(select(data.predictA1T12, Xnames))
  m.X.A0 <- apply(XmatA0, 2, mean)
  m.X.A1 <- apply(XmatA1, 2, mean)
  m.X.A0T12 <- apply(XmatA0T12, 2, mean)
  m.X.A1T12 <- apply(XmatA1T12, 2, mean)
  Xmat <- as.matrix(select(data, Xnames))
  # Calculate exp(X\beta) only once.
  exp.b001 <- exp(Xmat %*% coef(res$fit.list$fit.a0.01))
  exp.b002 <- exp(Xmat %*% coef(res$fit.list$fit.a0.02))
  exp.b012 <- exp(Xmat %*% coef(res$fit.list$fit.a0.12))
  exp.b101 <- exp(Xmat %*% coef(res$fit.list$fit.a1.01))
  exp.b102 <- exp(Xmat %*% coef(res$fit.list$fit.a1.02))
  exp.b112 <- exp(Xmat %*% coef(res$fit.list$fit.a1.12))
  # Sample bivariate fraily from gamma distribution
  theta.est <- (res$thetas[1] * sum(data$A==0) +
                  res$thetas[2] * sum(data$A==1))/nrow(data)
  gamma.scale <- theta.est
  gamma.common.shape <- rho/theta.est
  gamma.each.shape <-  (1 - rho)/theta.est
  gamma.common <- rgamma(n.sample.sim, shape = gamma.common.shape, scale = gamma.scale)
  gamma0 <- gamma.common + rgamma(n.sample.sim, shape = gamma.each.shape,
                                  scale = gamma.scale)
  gamma1 <- gamma.common + rgamma(n.sample.sim, shape = gamma.each.shape,
                                  scale = gamma.scale)
  # Sample X values from observed data
  ind.X <- sample(x = 1:n.sample, size = n.sample.sim, replace = T)
  exp.b001.sim <- exp.b001[ind.X, ]
  exp.b002.sim <- exp.b002[ind.X, ]
  exp.b012.sim <- exp.b012[ind.X, ]
  exp.b101.sim <- exp.b101[ind.X, ]
  exp.b102.sim <- exp.b102[ind.X, ]
  exp.b112.sim <- exp.b112[ind.X, ]
  step.A0T1 <- res$H.step.funcs$step.A0T1
  step.A0T2 <- res$H.step.funcs$step.A0T2
  step.A0T12 <- res$H.step.funcs$step.A0T12
  step.A1T1 <- res$H.step.funcs$step.A1T1
  step.A1T2 <- res$H.step.funcs$step.A1T2
  step.A1T12 <- res$H.step.funcs$step.A1T12
  # Replace values larger with tau with tau (later also fix probs)
  A0T1.times <- pmin(knots(step.A0T1), tau) %>% unique
  A0T2.times <- pmin(knots(step.A0T2), tau) %>% unique
  A0T12.times <- pmin(knots(step.A0T12), tau) %>% unique
  A1T1.times <- pmin(knots(step.A1T1), tau) %>% unique
  A1T2.times <- pmin(knots(step.A1T2), tau) %>% unique
  A1T12.times <- pmin(knots(step.A1T12), tau) %>% unique
  T1.0.sim <- T1.1.sim <- T2.0.sim <- T2.1.sim <- vector(length = n.sample.sim)
  # For the future: might be better to have linear interpolation between points using approxfun (but f gives global mehtod=constant).
  # For now, only sample from event times
  for(i in 1:n.sample.sim)
  {
  step.S.a0.01i <- function(t) {exp(-step.A0T1(t) * exp.b001.sim[i] * gamma0[i])}
  step.S.a0.02i <- function(t) {exp(-step.A0T2(t) * exp.b002.sim[i] * gamma0[i])}
  step.S.a0.12i <- function(t) {exp(-step.A0T12(t) * exp.b012.sim[i] * gamma0[i])}
  #H.a0.12.T1 <- step.A0T12(data.predictA0T12$T1) * exp.b012.sim[i]
  step.S.a1.01i <- function(t) {exp(-step.A1T1(t) * exp.b101.sim[i] * gamma1[i])}
  step.S.a1.02i <- function(t) {exp(-step.A1T2(t) * exp.b102.sim[i] * gamma1[i])}
  step.S.a1.12i <- function(t) {exp(-step.A1T12(t) * exp.b112.sim[i] * gamma1[i])}
  #H.a1.12.T1 <- step.A1T12(data.predictA1T12$T1) * exp.b112.sim[i]
  # Assign remaining prob mass to tau
  pr.A0T1 <- -diff(c(1, step.S.a0.01i(A0T1.times)))
  pr.A0T1[length(pr.A0T1)] <- pr.A0T1[length(pr.A0T1)] + 1 - sum(-diff(c(1, step.S.a0.01i(A0T1.times))))
  pr.A0T2 <- -diff(c(1, step.S.a0.02i(A0T2.times)))
  pr.A0T2[length(pr.A0T2)] <- pr.A0T2[length(pr.A0T2)] + 1 - sum(-diff(c(1, step.S.a0.02i(A0T2.times))))
  pr.A1T1 <- -diff(c(1, step.S.a1.01i(A1T1.times)))
  pr.A1T1[length(pr.A1T1)] <- pr.A1T1[length(pr.A1T1)] + 1 - sum(-diff(c(1, step.S.a1.01i(A1T1.times))))
  pr.A1T2 <- -diff(c(1, step.S.a1.02i(A1T2.times)))
  pr.A1T2[length(pr.A1T2)] <- pr.A1T2[length(pr.A1T2)] + 1 - sum(-diff(c(1, step.S.a1.02i(A1T2.times))))
  #### Simulate T1, T2

  T1.0.sim[i] <- sample(A0T1.times, 1, replace = T, prob = pr.A0T1)
  T2.0.sim[i] <- sample(A0T2.times, 1, replace = T, prob = pr.A0T2)
  T1.1.sim[i] <- sample(A1T1.times, 1, replace = T, prob = pr.A1T1)
  T2.1.sim[i] <- sample(A1T2.times, 1, replace = T, prob = pr.A1T2)
  if (T1.0.sim[i] <= T2.0.sim[i])
  {
    step.S.a0.12i <- function(t) {exp(-(step.A0T12(t) - step.A0T12(T1.0.sim[i])) * exp.b012.sim[i] * gamma0[i])}
    A0T12.times.i <- A0T12.times[A0T12.times > T1.0.sim[i]]
    pr.A0T12 <- -diff(c(1, step.S.a0.12i(A0T12.times.i)))
    pr.A0T12[length(pr.A0T12)] <- pr.A0T12[length(pr.A0T12)] + 1 - sum(-diff(c(1, step.S.a0.12i(A0T12.times.i))))
    T2.0.sim[i] <- sample(A0T12.times.i, 1, replace = T, prob = pr.A0T12)
  } else {
    T1.0.sim[i] <- Inf
  }
  if (T1.1.sim[i] <= T2.1.sim[i])
  {
    step.S.a1.12i <- function(t) {exp(-(step.A1T12(t) - step.A1T12(T1.1.sim[i])) * exp.b112.sim[i] * gamma1[i])}
    A1T12.times.i <- A1T12.times[A1T12.times > T1.1.sim[i]]
    pr.A1T12 <- -diff(c(1, step.S.a1.12i(A1T12.times.i)))
    pr.A1T12[length(pr.A1T12)] <- pr.A1T12[length(pr.A1T12)] + 1 - sum(-diff(c(1, step.S.a1.12i(A1T12.times.i))))
    T2.1.sim[i] <- sample(A1T12.times.i, 1, replace = T, prob = pr.A1T12)
  } else {
    T1.1.sim[i] <- Inf
  }
  }
  #### Estimated causal effects
  ad <- (T1.0.sim < Inf) & (T1.1.sim < Inf)
  nd <- (T1.0.sim == Inf) & (T1.1.sim == Inf)
    ########################################################
  ATE.T2.ad <- mean(T2.1.sim[ad] - T2.0.sim[ad])
  ATE.T2.nd <- mean(T2.1.sim[nd] - T2.0.sim[nd])
  ATE.T1.ad <- mean(T1.1.sim[ad] - T1.0.sim[ad])
  med.ATE.T2.ad <- median(T2.1.sim[ad]) - median(T2.0.sim[ad])
  med.ATE.T2.nd <- median(T2.1.sim[nd]) - median(T2.0.sim[nd])
  med.ATE.T1.ad <- median(T1.1.sim[ad]) - median(T1.0.sim[ad])
  list.ret <- list(ATE.T2.ad = ATE.T2.ad, ATE.T2.nd = ATE.T2.nd, ATE.T1.ad = ATE.T1.ad,
                   med.ATE.T2.ad = ATE.T2.ad, med.ATE.T2.nd = ATE.T2.nd, med.ATE.T1.ad = med.ATE.T1.ad)
  return(list.ret)
}
