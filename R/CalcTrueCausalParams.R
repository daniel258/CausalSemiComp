###########################################################################
# CasualSemiComp

# A function to calculate true parameter values
#Data-generating functions for simulations
###########################################################################
CalcTrueCausalParams <- function(n.sample, params, all.times = NULL, no.large, no.protected, adjusted = F,
                                 tau = NULL, X = NULL, RMST.only = F)
{
if(is.null(tau) & RMST.only==T) {stop("Only RMST is requested, but no tau is supplied")}
if(n.sample < 10001) {warning(paste0("Hard to believe that n.sample = ", n.sample, "is sufficient to
                                     approximate infinite population"))}
if (is.null(X))
  {
  sim.df <- SimDataWeibFrail(n.sample = n.sample, params, no.large = no.large,
                           no.protected = no.protected, cens.exp.rate = 0.000001)
} else {
  sim.df <- SimDataWeibFrail(n.sample = n.sample, params, no.large = no.large, X = X,
                             no.protected = no.protected, cens.exp.rate = 0.000001)
  }
T1.0 <- sim.df$T1.0
T1.1 <- sim.df$T1.1
T2.0 <- sim.df$T2.0
T2.1 <- sim.df$T2.1
true.prop.ad <- mean(T1.0 <= T2.0 & T1.1 <= T2.1)
true.prop.nd <- mean(T1.0 > T2.0 &  T1.1 > T2.1)
true.prop.dh <- mean(T1.0 > T2.0 &  T1.1 < T2.1)
if(RMST.only == F)
{
true.eta0 <- mean(T1.0 <= T2.0)
true.eta1 <- mean(T1.1 <= T2.1)


############################################################################################################
e0 <- ecdf(T2.0)
e1 <- ecdf(T2.1)
true.S2A0 <- 1 - e0(all.times)
true.S2A1 <- 1 - e1(all.times)
##### \eta_{A=a, T_2 \le t}
true.etas.T2leT.0 <- sapply(all.times, function(t) {mean(T1.0[T2.0 <= t] <= T2.0[T2.0 <= t])})
true.etas.T2leT.1 <- sapply(all.times, function(t) {mean(T1.1[T2.1 <= t] <= T2.1[T2.1 <= t])})

T1.0.inf <- T1.0
T1.0.inf[T1.0 > T2.0] <- Inf
T1.1.inf <- T1.1
T1.1.inf[T1.1 > T2.1] <- Inf
true.S1A0 <- sapply(all.times, function(t) {mean(T1.0.inf > t)})
true.S1A1 <- sapply(all.times, function(t) {mean(T1.1.inf > t)})
##### S_{1|A=a, T_1 \le T_2}
true.S1A0T1leT2 <- sapply(all.times, function(t) {mean(T1.0[T1.0 < T2.0] > t)})
true.S1A1T1leT2 <- sapply(all.times, function(t) {mean(T1.1[T1.1 < T2.1] > t)})
##### S_{2|A=a, T_1 \le T_2}
true.S2A0T1leT2 <- sapply(all.times, function(t) {mean(T2.0[T1.0 < T2.0] > t)})
true.S2A1T1leT2 <- sapply(all.times, function(t) {mean(T2.1[T1.1 < T2.1] > t)})
#####  S_{2|A=a, T_1 > T_2}
true.S2A0T1gT2 <- sapply(all.times, function(t) {mean(T2.0[T1.0 > T2.0] > t)})
true.S2A1T1gT2 <- sapply(all.times, function(t) {mean(T2.1[T1.1 > T2.1] > t)})
########################################################
F2.a0.ad <- sapply(all.times, function(t) {mean(T2.0[T1.0 < T2.0 & T1.1 < T2.1] <= t)})
F2.a1.ad <- sapply(all.times, function(t) {mean(T2.1[T1.0 < T2.0 & T1.1 < T2.1] <= t)})
F2.a0.nd <- sapply(all.times, function(t) {mean(T2.0[T1.0 > T2.0 & T1.1 > T2.1] <= t)})
F2.a1.nd <- sapply(all.times, function(t) {mean(T2.1[T1.0 > T2.0 & T1.1 > T2.1] <= t)})
F1.a0.ad <- sapply(all.times, function(t) {mean(T1.0.inf[T1.0 < T2.0 & T1.1 < T2.1] <= t)})
F1.a1.ad <- sapply(all.times, function(t) {mean(T1.1.inf[T1.0 < T2.0 & T1.1 < T2.1] <= t)})
F1.a1.h <- sapply(all.times, function(t) {mean(T1.1.inf[T1.0 > T2.0 & T1.1 < T2.1] <= t)})
if (adjusted==T)
{
  Z <- sim.df$X[, 1]
  Zprob.ad <- mean(Z[T1.0 < T2.0])
  Zprob.nd <- mean(Z[T1.1 > T2.1])
  true.eta0Z0 <- mean(T1.0[Z==0] <= T2.0[Z==0])
  true.eta1Z0 <- mean(T1.1[Z==0] <= T2.1[Z==0])
  true.eta0Z1 <- mean(T1.0[Z==1] <= T2.0[Z==1])
  true.eta1Z1 <- mean(T1.1[Z==1] <= T2.1[Z==1])

  ############################################################################################################
  e0Z0 <- ecdf(T2.0[Z==0])
  e1Z0 <- ecdf(T2.1[Z==0])
  e0Z1 <- ecdf(T2.0[Z==1])
  e1Z1 <- ecdf(T2.1[Z==1])
  true.S2A0Z0 <- 1 - e0Z0(all.times)
  true.S2A1Z0 <- 1 - e1Z0(all.times)
  true.S2A0Z1 <- 1 - e0Z1(all.times)
  true.S2A1Z1 <- 1 - e1Z1(all.times)
  ##### \eta_{A=a, T_2 \le t}
  true.etas.T2leT.0.Z0 <- sapply(all.times, function(t) {mean(T1.0[T2.0 <= t & Z==0] <=
                                                               T2.0[T2.0 <= t & Z==0])})
  true.etas.T2leT.1.Z0 <- sapply(all.times, function(t) {mean(T1.1[T2.1 <= t & Z==0] <=
                                                               T2.1[T2.1 <= t & Z==0])})
  true.etas.T2leT.0.Z1 <- sapply(all.times, function(t) {mean(T1.0[T2.0 <= t & Z==1] <=
                                                                T2.0[T2.0 <= t & Z==1])})
  true.etas.T2leT.1.Z1 <- sapply(all.times, function(t) {mean(T1.1[T2.1 <= t & Z==1] <=
                                                                T2.1[T2.1 <= t & Z==1])})

  T1.0.inf <- T1.0
  T1.0.inf[T1.0 > T2.0] <- Inf
  T1.1.inf <- T1.1
  T1.1.inf[T1.1 > T2.1] <- Inf
  true.S1A0Z0 <- sapply(all.times, function(t) {mean(T1.0.inf[Z==0] > t)})
  true.S1A1Z0 <- sapply(all.times, function(t) {mean(T1.1.inf[Z==0] > t)})
  true.S1A0Z1 <- sapply(all.times, function(t) {mean(T1.0.inf[Z==1] > t)})
  true.S1A1Z1 <- sapply(all.times, function(t) {mean(T1.1.inf[Z==1] > t)})
  ##### S_{1|A=a, T_1 \le T_2}
  true.S1A0T1leT2Z0 <- sapply(all.times, function(t) {mean(T1.0[T1.0 < T2.0 & Z==0] > t)})
  true.S1A1T1leT2Z0 <- sapply(all.times, function(t) {mean(T1.1[T1.1 < T2.1 & Z==0] > t)})
  true.S1A0T1leT2Z1 <- sapply(all.times, function(t) {mean(T1.0[T1.0 < T2.0 & Z==1] > t)})
  true.S1A1T1leT2Z1 <- sapply(all.times, function(t) {mean(T1.1[T1.1 < T2.1 & Z==1] > t)})
  ##### S_{2|A=a, T_1 \le T_2}
  true.S2A0T1leT2Z0 <- sapply(all.times, function(t) {mean(T2.0[T1.0 < T2.0 & Z==0] > t)})
  true.S2A1T1leT2Z0 <- sapply(all.times, function(t) {mean(T2.1[T1.1 < T2.1 & Z==0] > t)})
  true.S2A0T1leT2Z1 <- sapply(all.times, function(t) {mean(T2.0[T1.0 < T2.0 & Z==1] > t)})
  true.S2A1T1leT2Z1 <- sapply(all.times, function(t) {mean(T2.1[T1.1 < T2.1 & Z==1] > t)})
  #####  S_{2|A=a, T_1 > T_2}
  true.S2A0T1gT2Z0 <- sapply(all.times, function(t) {mean(T2.0[T1.0 > T2.0 & Z==0] > t)})
  true.S2A1T1gT2Z0 <- sapply(all.times, function(t) {mean(T2.1[T1.1 > T2.1 & Z==0] > t)})
  true.S2A0T1gT2Z1 <- sapply(all.times, function(t) {mean(T2.0[T1.0 > T2.0 & Z==1] > t)})
  true.S2A1T1gT2Z1 <- sapply(all.times, function(t) {mean(T2.1[T1.1 > T2.1 & Z==1] > t)})
  ######
  F2.a0.Z0.ad <- sapply(all.times, function(t) {mean(T2.0[T1.0 < T2.0 & T1.1 < T2.1 & Z==0] <= t)})
  F2.a1.Z0.ad <- sapply(all.times, function(t) {mean(T2.1[T1.0 < T2.0 & T1.1 < T2.1 & Z==0] <= t)})
  F2.a0.Z1.ad <- sapply(all.times, function(t) {mean(T2.0[T1.0 < T2.0 & T1.1 < T2.1 & Z==1] <= t)})
  F2.a1.Z1.ad <- sapply(all.times, function(t) {mean(T2.1[T1.0 < T2.0 & T1.1 < T2.1 & Z==1] <= t)})
  F2.a0.Z0.nd <- sapply(all.times, function(t) {mean(T2.0[T1.0 > T2.0 & T1.1 > T2.1 & Z==0] <= t)})
  F2.a1.Z0.nd <- sapply(all.times, function(t) {mean(T2.1[T1.0 > T2.0 & T1.1 > T2.1 & Z==0] <= t)})
  F2.a0.Z1.nd <- sapply(all.times, function(t) {mean(T2.0[T1.0 > T2.0 & T1.1 > T2.1 & Z==1] <= t)})
  F2.a1.Z1.nd <- sapply(all.times, function(t) {mean(T2.1[T1.0 > T2.0 & T1.1 > T2.1 & Z==1] <= t)})
  F1.a0.Z0.ad <- sapply(all.times, function(t) {mean(T1.0.inf[T1.0 < T2.0 & T1.1 < T2.1 & Z==0] <= t)})
  F1.a1.Z0.ad <- sapply(all.times, function(t) {mean(T1.1.inf[T1.0 < T2.0 & T1.1 < T2.1 & Z==0] <= t)})
  F1.a0.Z1.ad <- sapply(all.times, function(t) {mean(T1.0.inf[T1.0 < T2.0 & T1.1 < T2.1 & Z==1] <= t)})
  F1.a1.Z1.ad <- sapply(all.times, function(t) {mean(T1.1.inf[T1.0 < T2.0 & T1.1 < T2.1 & Z==1] <= t)})
  F1.a1.Z0.h <- sapply(all.times, function(t) {mean(T1.1.inf[T1.0 > T2.0 & T1.1 < T2.1 & Z==0] <= t)})
  F1.a1.Z1.h <- sapply(all.times, function(t) {mean(T1.1.inf[T1.0 > T2.0 & T1.1 < T2.1 & Z==1] <= t)})
  #### Create list to save everything #####
  true.eta0 <- list(true.eta0 = true.eta0, true.eta0Z0 = true.eta0Z0, true.eta0Z1 = true.eta0Z1)
  true.eta1 <- list(true.eta1 = true.eta1, true.eta1Z0 = true.eta1Z0, true.eta1Z1 = true.eta1Z1)
  true.S2A0 <- list(true.S2A0 = true.S2A0, true.S2A0Z0 = true.S2A0Z0, true.S2A0Z1 = true.S2A0Z1)
  true.S2A1 <- list(true.S2A1 = true.S2A1, true.S2A1Z0 = true.S2A1Z0, true.S2A1Z1 = true.S2A1Z1)
  true.S1A0 <- list(true.S1A0 = true.S1A0, true.S1A0Z0 = true.S1A0Z0, true.S1A0Z1 = true.S1A0Z1)
  true.S1A1 <- list(true.S1A1 = true.S1A1, true.S1A1Z0 = true.S1A1Z0, true.S1A1Z1 = true.S1A1Z1)
  true.etas.T2leT.0 <- list(true.etas.T2leT.0 = true.etas.T2leT.0,
                            true.etas.T2leT.0.Z0 = true.etas.T2leT.0.Z0,
                            true.etas.T2leT.0.Z1 = true.etas.T2leT.0.Z1)
  true.etas.T2leT.1 <- list(true.etas.T2leT.1 = true.etas.T2leT.1,
                            true.etas.T2leT.1.Z0 = true.etas.T2leT.1.Z0,
                            true.etas.T2leT.1.Z1 = true.etas.T2leT.1.Z1)
  true.S1A0T1leT2 <- list(true.S1A0T1leT2 = true.S1A0T1leT2, true.S1A0T1leT2Z0 = true.S1A0T1leT2Z0,
                          true.S1A0T1leT2Z1 = true.S1A0T1leT2Z1)
  true.S1A1T1leT2 <- list(true.S1A1T1leT2 = true.S1A1T1leT2, true.S1A1T1leT2Z0 = true.S1A1T1leT2Z0,
                          true.S1A1T1leT2Z1 = true.S1A1T1leT2Z1)
  true.S2A0T1leT2 <- list(true.S2A0T1leT2 = true.S2A0T1leT2, true.S2A0T1leT2Z0 = true.S2A0T1leT2Z0,
                          true.S2A0T1leT2Z1 = true.S2A0T1leT2Z1)
  true.S2A1T1leT2 <- list(true.S2A1T1leT2 = true.S2A1T1leT2, true.S2A1T1leT2Z0 = true.S2A1T1leT2Z0,
                          true.S2A1T1leT2Z1 = true.S2A1T1leT2Z1)
  true.S2A0T1gT2 <- list(true.S2A0T1gT2 = true.S2A0T1gT2, true.S2A0T1gT2Z0 = true.S2A0T1gT2Z0,
                         true.S2A0T1gT2Z1 = true.S2A0T1gT2Z1)
  true.S2A1T1gT2 <- list(true.S2A1T1gT2 = true.S2A1T1gT2, true.S2A1T1gT2Z0 = true.S2A1T1gT2Z0,
                         true.S2A1T1gT2Z1 = true.S2A1T1gT2Z1)
  F2.a0.ad <- list(F2.a0.ad = F2.a0.ad, F2.a0.Z0.ad = F2.a0.Z0.ad, F2.a0.Z1.ad = F2.a0.Z1.ad)
  F2.a1.ad <- list(F2.a1.ad = F2.a1.ad, F2.a1.Z0.ad = F2.a1.Z0.ad, F2.a1.Z1.ad = F2.a1.Z1.ad)
  F2.a0.nd <- list(F2.a0.nd = F2.a0.nd, F2.a0.Z0.nd = F2.a0.Z0.nd, F2.a0.Z1.nd = F2.a0.Z1.nd)
  F2.a1.nd <- list(F2.a1.nd = F2.a1.nd, F2.a1.Z0.nd = F2.a1.Z0.nd, F2.a1.Z1.nd = F2.a1.Z1.nd)
  F1.a0.ad <- list(F1.a0.ad = F1.a0.ad, F1.a0.Z0.ad = F1.a0.Z0.ad, F1.a0.Z1.ad = F1.a0.Z1.ad)
  F1.a1.ad <- list(F1.a1.ad = F1.a1.ad, F1.a1.Z0.ad = F1.a1.Z0.ad, F1.a1.Z1.ad = F1.a1.Z1.ad)
  F1.a1.h <- list(F1.a1.h = F1.a1.h, F1.a1.Z0.h = F1.a1.Z0.h, F1.a1.Z1.h = F1.a1.Z1.h)

list.ret <- list(true.eta0 = true.eta0, true.eta1 = true.eta1,
                 true.prop.ad = true.prop.ad, true.prop.nd = true.prop.nd, true.prop.dh = true.prop.dh,
                 true.S2A0 = true.S2A0, true.S2A1 = true.S2A1,
                 true.S1A0 = true.S1A0, true.S1A1 = true.S1A1,
                 true.etas.T2leT.0 = true.etas.T2leT.0, true.etas.T2leT.1 = true.etas.T2leT.1,
                 true.S1A0T1leT2 = true.S1A0T1leT2, true.S1A1T1leT2 = true.S1A1T1leT2,
                 true.S2A0T1leT2 = true.S2A0T1leT2, true.S2A1T1leT2 = true.S2A1T1leT2,
                 true.S2A0T1gT2 = true.S2A0T1gT2, true.S2A1T1gT2 = true.S2A1T1gT2,
                 F2.a0.ad = F2.a0.ad, F2.a1.ad = F2.a1.ad, F2.a0.nd = F2.a0.nd, F2.a1.nd = F2.a1.nd,
                 F1.a0.ad = F1.a0.ad, F1.a1.ad = F1.a1.ad, F1.a1.h = F1.a1.h,
                 Zprob.ad = Zprob.ad, Zprob.nd = Zprob.nd,
                 ATE.T2.ad = NULL, ATE.T2.nd = NULL, ATE.T1.ad = NULL,
                 med.ATE.T2.ad = NULL, med.ATE.T2.nd = NULL, med.ATE.T1.ad = NULL)
} else {
  list.ret <- list(true.eta0 = true.eta0, true.eta1 = true.eta1,
                   true.prop.ad = true.prop.ad, true.prop.nd = true.prop.nd, true.prop.dh = true.prop.dh,
                   true.S2A0 = true.S2A0, true.S2A1 = true.S2A1,
                   true.S1A0 = true.S1A0, true.S1A1 = true.S1A1,
                   true.etas.T2leT.0 = true.etas.T2leT.0, true.etas.T2leT.1 = true.etas.T2leT.1,
                   true.S1A0T1leT2 = true.S1A0T1leT2, true.S1A1T1leT2 = true.S1A1T1leT2,
                   true.S2A0T1leT2 = true.S2A0T1leT2, true.S2A1T1leT2 = true.S2A1T1leT2,
                   true.S2A0T1gT2 = true.S2A0T1gT2, true.S2A1T1gT2 = true.S2A1T1gT2,
                   F2.a0.ad = F2.a0.ad, F2.a1.ad = F2.a1.ad, F2.a0.nd = F2.a0.nd, F2.a1.nd = F2.a1.nd,
                   F1.a0.ad = F1.a0.ad, F1.a1.ad = F1.a1.ad, F1.a1.h = F1.a1.h,
                   ATE.T2.ad = NULL, ATE.T2.nd = NULL, ATE.T1.ad = NULL,
                   med.ATE.T2.ad = NULL, med.ATE.T2.nd = NULL, med.ATE.T1.ad = NULL)
}}
else {
  #if RMST.only==T, all nonparam effects are not computed and are not computed nor returned
  list.ret <- list(true.prop.ad = true.prop.ad, true.prop.nd = true.prop.nd, true.prop.dh = true.prop.dh)
}
# add RMST if tau is non null
if (!is.null(tau))
{
  T1.0.tau <- pmin(T1.0, tau)
  T1.1.tau <- pmin(T1.1, tau)
  T2.0.tau <- pmin(T2.0, tau)
  T2.1.tau <- pmin(T2.1, tau)
  ad <- (T1.0.tau < T2.0.tau) & (T1.1.tau < T2.1.tau)
  nd <- (T1.0.tau >= T2.0.tau) & (T1.1.tau >= T2.1.tau)
  prop.ad.tau <- mean(ad)
  prop.nd.tau <- mean(nd)
  ########################################################
  mean.T2.ad.a1 <- mean(T2.1.tau[ad])
  mean.T2.ad.a0 <- mean(T2.0.tau[ad])
  mean.T1.ad.a1 <- mean(T1.1.tau[ad])
  mean.T1.ad.a0 <- mean(T1.0.tau[ad])
  mean.T2.nd.a1 <- mean(T2.1.tau[nd])
  mean.T2.nd.a0 <- mean(T2.0.tau[nd])
  med.T2.ad.a1 <- median(T2.1.tau[ad])
  med.T2.ad.a0 <- median(T2.0.tau[ad])
  med.T1.ad.a1 <- median(T1.1.tau[ad])
  med.T1.ad.a0 <- median(T1.0.tau[ad])
  med.T2.nd.a1 <- median(T2.1.tau[nd])
  med.T2.nd.a0 <- median(T2.0.tau[nd])
  ATE.T2.ad <- mean(T2.1.tau[ad] - T2.0.tau[ad])
  ATE.T2.nd <- mean(T2.1.tau[nd] - T2.0.tau[nd])
  ATE.T1.ad <- mean(T1.1.tau[ad] - T1.0.tau[ad])
  med.ATE.T2.ad <- median(T2.1.tau[ad]) - median(T2.0.tau[ad])
  med.ATE.T2.nd <- median(T2.1.tau[nd]) - median(T2.0.tau[nd])
  med.ATE.T1.ad <- median(T1.1.tau[ad]) - median(T1.0.tau[ad])
  list.ret$prop.ad.tau <- prop.ad.tau
  list.ret$prop.nd.tau <- prop.nd.tau
  list.ret$mean.T2.ad.a1 <- mean.T2.ad.a1
  list.ret$mean.T2.ad.a0 <- mean.T2.ad.a0
  list.ret$mean.T1.ad.a1 <- mean.T1.ad.a1
  list.ret$mean.T1.ad.a0 <- mean.T1.ad.a0
  list.ret$mean.T2.nd.a1 <- mean.T2.nd.a1
  list.ret$mean.T2.nd.a0 <- mean.T2.nd.a0
  list.ret$med.T2.ad.a1 <- med.T2.ad.a1
  list.ret$med.T2.ad.a0 <- med.T2.ad.a0
  list.ret$med.T1.ad.a1 <- med.T1.ad.a1
  list.ret$med.T1.ad.a0 <- med.T1.ad.a0
  list.ret$med.T2.nd.a1 <- med.T2.nd.a1
  list.ret$med.T2.nd.a0 <- med.T2.nd.a0
  list.ret$ATE.T2.ad <- ATE.T2.ad
  list.ret$ATE.T2.nd <- ATE.T2.nd
  list.ret$ATE.T1.ad <- ATE.T1.ad
  list.ret$med.ATE.T2.ad <- med.ATE.T2.ad
  list.ret$med.ATE.T2.nd <- med.ATE.T2.nd
  list.ret$med.ATE.T1.ad <- med.ATE.T1.ad
  # list.ret$T1.0 = T1.0; list.ret$T2.0 = T2.0
  # list.ret$T1.1 = T1.1; list.ret$T2.1 = T2.1
  # list.ret$T1.0.tau = T1.0.tau; list.ret$T2.0.tau = T2.0.tau
  # list.ret$T1.1.tau = T1.1.tau; list.ret$T2.1.tau = T2.1.tau
}
########################################################
return(list.ret)
}
