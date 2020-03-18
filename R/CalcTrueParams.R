###########################################################################
# CasualSemiComp

# A function to calculate true parameter values
#Data-generating functions for simulations
###########################################################################
CalcTrueParams <- function(n.sample,
                         base.weib.scale.a0.01, base.weib.scale.a1.01,
                         base.weib.scale.a0.02, base.weib.scale.a1.02,
                         base.weib.scale.a0.12, base.weib.scale.a1.12,
                         base.weib.shape.a0.01, base.weib.shape.a1.01,
                         base.weib.shape.a0.02, base.weib.shape.a1.02,
                         base.weib.shape.a0.12, base.weib.shape.a1.12,
                         theta,
                         beta.a0.01, beta.a1.01,
                         beta.a0.02, beta.a1.02,
                         beta.a0.12, beta.a1.12,
                         all.times)
{
if(n.sample < 10001) {warning(paste0("Hard to believe that n.sample = ", n.sample, "is sufficient to
                                     approximate infinite population"))}
sim.df <- SimDataWeibFrailNoDP(n.sample = n.sample,
                               base.weib.shape.a0.01 = base.weib.shape.a0.01,
                               base.weib.shape.a0.02 = base.weib.shape.a0.02,
                               base.weib.shape.a0.12 = base.weib.shape.a0.12,
                               base.weib.shape.a1.01 = base.weib.shape.a1.01,
                               base.weib.shape.a1.02 = base.weib.shape.a1.02,
                               base.weib.shape.a1.12 = base.weib.shape.a1.12,
                               base.weib.scale.a0.01 = base.weib.scale.a0.01,
                               base.weib.scale.a0.02 = base.weib.scale.a0.02,
                               base.weib.scale.a0.12 = base.weib.scale.a0.12,
                               base.weib.scale.a1.01 = base.weib.scale.a1.01,
                               base.weib.scale.a1.02 = base.weib.scale.a1.02,
                               base.weib.scale.a1.12 = base.weib.scale.a1.12,
                               theta = theta,
                               beta.a0.01 = beta.a0.01, beta.a1.01 = beta.a1.01,
                               beta.a0.02 = beta.a0.02, beta.a1.02 = beta.a1.02,
                               beta.a0.12 = beta.a0.12, beta.a1.12 = beta.a1.12, cens.exp.rate = 0.000001)
AD.0 <- sim.df$AD.0
AD.1 <- sim.df$AD.1
Death.0 <- sim.df$Death.0
Death.1 <- sim.df$Death.1

true.eta0 <- mean(AD.0 < Death.0)
true.eta1 <- mean(AD.1 < Death.1)
true.prop.ad <- mean(AD.0 <= Death.0 & AD.1 <= Death.1)
true.prop.nd <- mean(AD.0 > Death.0 &  AD.1 > Death.1)
true.prop.dh <- 1 - true.prop.ad - true.prop.nd

# true.eta0 <- mean(AD.0 <= Death.0 & Death.0 <= 100)
# true.eta1 <- mean(AD.1 <= Death.1 & Death.1 <= 100)
# true.prop.ad <- mean(AD.0 <= Death.0 & AD.1 <= Death.1 & Death.0 <= 100 & Death.1 <= 100)
# true.prop.nd <- mean(AD.0 > pmin(Death.0, 100) & pmin(AD.1 > Death.1, 100))
# true.prop.dh <- 1 - true.prop.ad - true.prop.nd

############################################################################################################
e0 <- ecdf(Death.0)
e1 <- ecdf(Death.1)
true.S2A0 <- 1 - e0(all.times)
true.S2A1 <- 1 - e1(all.times)
##### \eta_{A=a, T_2 \le t}
# true.etas.T2leT.0 <- sapply(all.times, function(t) {mean(AD.0[Death.0 <= t] <= Death.0[Death.0 <= t])})*mean(Death.0 <= 100)
# true.etas.T2leT.1 <- sapply(all.times, function(t) {mean(AD.1[Death.1 <= t] <= Death.1[Death.1 <= t])})*mean(Death.1 <= 100)
true.etas.T2leT.0 <- sapply(all.times, function(t) {mean(AD.0[Death.0 <= t] <= Death.0[Death.0 <= t])})
true.etas.T2leT.1 <- sapply(all.times, function(t) {mean(AD.1[Death.1 <= t] <= Death.1[Death.1 <= t])})

AD.0.inf <- AD.0
AD.0.inf[AD.0 > Death.0] <- Inf
AD.1.inf <- AD.1
AD.1.inf[AD.1 > Death.1] <- Inf
true.S1A0 <- sapply(all.times, function(t) {mean(AD.0.inf > t)})
true.S1A1 <- sapply(all.times, function(t) {mean(AD.1.inf > t)})
##### S_{1|A=a, T_1 \le T_2}
true.S1A0T1leT2 <- sapply(all.times, function(t) {mean(AD.0[AD.0 < Death.0] > t)})
true.S1A1T1leT2 <- sapply(all.times, function(t) {mean(AD.1[AD.1 < Death.1] > t)})
##### S_{2|A=a, T_1 \le T_2}
true.S2A0T1leT2 <- sapply(all.times, function(t) {mean(Death.0[AD.0 < Death.0] > t)})
true.S2A1T1leT2 <- sapply(all.times, function(t) {mean(Death.1[AD.1 < Death.1] > t)})
#####  S_{2|A=a, T_1 > T_2}
true.S2A0T1gT2 <- sapply(all.times, function(t) {mean(Death.0[AD.0 > Death.0] > t)})
true.S2A1T1gT2 <- sapply(all.times, function(t) {mean(Death.1[AD.1 > Death.1] > t)})
list.ret <- list(true.eta0 = true.eta0, true.eta1 = true.eta1,
                 true.prop.ad = true.prop.ad, true.prop.nd = true.prop.nd, true.prop.dh = true.prop.dh,
                 true.S2A0 = true.S2A0, true.S2A1 = true.S2A1,
                 true.S1A0 = true.S1A0, true.S1A1 = true.S1A1,
                 true.etas.T2leT.0 = true.etas.T2leT.0, true.etas.T2leT.1 = true.etas.T2leT.1,
                 true.S1A0T1leT2 = true.S1A0T1leT2, true.S1A1T1leT2 = true.S1A1T1leT2,
                 true.S2A0T1leT2 = true.S2A0T1leT2, true.S2A1T1leT2 = true.S2A1T1leT2,
                 true.S2A0T1gT2 = true.S2A0T1gT2, true.S2A1T1gT2 = true.S2A1T1gT2)
}
