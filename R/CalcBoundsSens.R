########################################################################
## Function for CausalSemiComp
### Estimation of various componenets for bounds and senstivity analysis

### The function returns the following:
#  S_{2|A=a}(t), for a=0,1  (names: S2A0 and S2A1)
#  \eta_{A=a}, for a=0,1  (names: etaA0 and etaA1)
#  \eta_{A=a, T_2 \le t} for a=0,1 (names: etasA0T2.le.t and etasA1T2.le.t)
#  S_{1|A=a}(t) for a=0,1 (names: S1A0.all.times and S1A1.all.times)
#  all.times
########################################################################
CalcBoundsSens <- function(surv.comp, bounds = T, sens = F, sens.param = NULL, probs = T, RMST = T)
{
  S2A0 <- surv.comp$S2A0
  S2A1 <- surv.comp$S2A1
  etaA0 <- surv.comp$etaA0
  etaA1 <- surv.comp$etaA1
  etasA0T2.le.t <- surv.comp$etasA0T2.le.t
  etasA1T2.le.t <- surv.comp$etasA1T2.le.t
  S1A0.all.times <- surv.comp$S1A0.all.times
  S1A1.all.times <- surv.comp$S1A1.all.times
  all.times <- surv.comp$all.times
  if(bounds==T)
  {
  F2A0T1lT2 <- (1 - surv.comp$S2A0)*etasA0T2.le.t/etaA0
  F1A0T1lT2 <- (1 - S1A0.all.times)*etasA0T2.le.t/etaA0
  ad.T2.L <- pmax(0, 1 - S2A1$surv/etaA0) - F2A0T1lT2
  ad.T2.H <- pmin(1, (1 - S2A1$surv) * etasA1T2.le.t/etaA0) - F2A0T1lT2
  nd.T2.L <- MISSING - pmin(1, (1 - S2A0$surv) * (1 - etasA0T2.le.t) / (1 - etaA1))
  nd.T2.H <- MISSING - pmax(0, 1 - (1 - S2A0$surv) / (1 - etaA1))
  #ad.T1.L <-
  #ad.T1.H <-
    ## bound for ad CDF

  }
}
