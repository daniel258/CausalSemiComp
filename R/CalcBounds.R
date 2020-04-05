########################################################################
## Function for CausalSemiComp
### Estimation of various componenets for bounds and senstivity analysis

### The function returns the following:
########################################################################
CalcBounds <- function(Res, Z = NULL)
{
  if (is.null(Z))
    {
    all.times <- Res$all.times
    S2A0.all.times <- Res$S2A0.all.times
    S2A1.all.times <- Res$S2A1.all.times
    etaA0 <- Res$etaA0
    etaA1 <- Res$etaA1
    etasA0T2.le.t <- Res$etasA0T2.le.t
    etasA1T2.le.t <- Res$etasA1T2.le.t
    S1A0.all.times <- Res$S1A0.all.times
    S1A1.all.times <- Res$S1A1.all.times
    S1A0T1lT2.all.times <- Res$S1A0T1lT2.all.times
    S1A1T1lT2.all.times <- Res$S1A1T1lT2.all.times
    F2A0T1lT2.all.times <- (1 - S2A0.all.times)*etasA0T2.le.t/etaA0
    F1A0T1lT2.all.times <- 1 - S1A0T1lT2.all.times
    F2A0T1gT2 <- (1 - etasA0T2.le.t)*(1 - S2A0.all.times)/(1 - etaA0)
    F2A1T1gT2 <- (1 - etasA1T2.le.t)*(1 - S2A1.all.times)/(1 - etaA1)
    ad.T2.L <- pmax(0, 1 - S2A1.all.times/etaA0) - F2A0T1lT2.all.times
    ad.T2.U <- pmin(1, (1 - S2A1.all.times) * etasA1T2.le.t/etaA0) - F2A0T1lT2.all.times
    nd.T2.L <- F2A1T1gT2 - pmin(1, (1 - S2A0.all.times) * (1 - etasA0T2.le.t) / (1 - etaA1))
    nd.T2.U <- F2A1T1gT2 - pmax(0, 1 - (1 - S2A0.all.times) / (1 - etaA1))
    ad.T1.L <- pmax(0, 1 - S1A1.all.times/etaA0) - F1A0T1lT2.all.times
    ad.T1.U <- pmin(1, (1 - S1A1.all.times)/etaA0) - F1A0T1lT2.all.times
    adjusted <- F
  } else {
    all.times <- Res$all.times
    ResZ0 <- Res$ResZ0
    ResZ1 <- Res$ResZ1
    Zprob.ad <- Res$Zprob.ad # This is the probability Pr(Z = 1| T_1(0)\le T_2(0))
    Zprob.nd <- Res$Zprob.nd # This is the probability Pr(Z = 1| T_1(1) > T_2(1))
    S2A0.all.timesZ0 <- ResZ0$S2A0.all.times
    S2A1.all.timesZ0 <- ResZ0$S2A1.all.times
    S2A0.all.timesZ1 <- ResZ1$S2A0.all.times
    S2A1.all.timesZ1 <- ResZ1$S2A1.all.times
    etaA0Z0 <- ResZ0$etaA0
    etaA1Z0 <- ResZ0$etaA1
    etaA0Z1 <- ResZ1$etaA0
    etaA1Z1 <- ResZ1$etaA1
    etasA0T2.le.tZ0 <- ResZ0$etasA0T2.le.t
    etasA1T2.le.tZ0 <- ResZ0$etasA1T2.le.t
    etasA0T2.le.tZ1 <- ResZ1$etasA0T2.le.t
    etasA1T2.le.tZ1 <- ResZ1$etasA1T2.le.t
    S1A0.all.timesZ0 <- ResZ0$S1A0.all.times
    S1A1.all.timesZ0 <- ResZ0$S1A1.all.times
    S1A0.all.timesZ1 <- ResZ1$S1A0.all.times
    S1A1.all.timesZ1 <- ResZ1$S1A1.all.times
    S1A0T1lT2.all.timesZ0 <- ResZ0$S1A0T1lT2.all.times
    S1A1T1lT2.all.timesZ0 <- ResZ0$S1A1T1lT2.all.times
    S1A0T1lT2.all.timesZ1 <- ResZ1$S1A0T1lT2.all.times
    S1A1T1lT2.all.timesZ1 <- ResZ1$S1A1T1lT2.all.times
    F2A0T1lT2.all.timesZ0 <- (1 - S2A0.all.timesZ0)*etasA0T2.le.tZ0/etaA0Z0
    F1A0T1lT2.all.timesZ0 <- 1 - S1A0T1lT2.all.timesZ0
    F2A0T1lT2.all.timesZ1 <- (1 - S2A0.all.timesZ1)*etasA0T2.le.tZ1/etaA0Z1
    F1A0T1lT2.all.timesZ1 <- 1 - S1A0T1lT2.all.timesZ1
    F2A0T1gT2Z0 <- (1 - etasA0T2.le.tZ0) * (1 - S2A0.all.timesZ0)/(1 - etaA0Z0)
    F2A1T1gT2Z0 <- (1 - etasA1T2.le.tZ0) * (1 - S2A1.all.timesZ0)/(1 - etaA1Z0)
    F2A0T1gT2Z1 <- (1 - etasA0T2.le.tZ1) * (1 - S2A0.all.timesZ1)/(1 - etaA0Z1)
    F2A1T1gT2Z1 <- (1 - etasA1T2.le.tZ1) * (1 - S2A1.all.timesZ1)/(1 - etaA1Z1)
    ad.T2.L.Z0 <- pmax(0, 1 - S2A1.all.timesZ0/etaA0Z0) - F2A0T1lT2.all.timesZ0
    ad.T2.U.Z0 <- pmin(1, (1 - S2A1.all.timesZ0) * etasA1T2.le.tZ0/etaA0Z0) - F2A0T1lT2.all.timesZ0
    ad.T2.L.Z1 <- pmax(0, 1 - S2A1.all.timesZ1/etaA0Z1) - F2A0T1lT2.all.timesZ1
    ad.T2.U.Z1 <- pmin(1, (1 - S2A1.all.timesZ1) * etasA1T2.le.tZ1/etaA0Z1) - F2A0T1lT2.all.timesZ1
    nd.T2.L.Z0 <- F2A1T1gT2Z0 - pmin(1, (1 - S2A0.all.timesZ0) * (1 - etasA0T2.le.tZ0) / (1 - etaA1Z0))
    nd.T2.U.Z0 <- F2A1T1gT2Z0 - pmax(0, 1 - (1 - S2A0.all.timesZ0) / (1 - etaA1Z0))
    nd.T2.L.Z1 <- F2A1T1gT2Z1 - pmin(1, (1 - S2A0.all.timesZ1) * (1 - etasA0T2.le.tZ1) / (1 - etaA1Z1))
    nd.T2.U.Z1 <- F2A1T1gT2Z1 - pmax(0, 1 - (1 - S2A0.all.timesZ1) / (1 - etaA1Z1))
    ad.T1.L.Z0 <- pmax(0, 1 - S1A1.all.timesZ0/etaA0Z0) - F1A0T1lT2.all.timesZ0
    ad.T1.L.Z1 <- pmax(0, 1 - S1A1.all.timesZ1/etaA0Z1) - F1A0T1lT2.all.timesZ1
    ad.T1.U.Z0 <- pmin(1, (1 - S1A1.all.timesZ0)/etaA0Z0) - F1A0T1lT2.all.timesZ0
    ad.T1.U.Z1 <- pmin(1, (1 - S1A1.all.timesZ1)/etaA0Z1) - F1A0T1lT2.all.timesZ1
    ad.T2.L <- (1 - Zprob.ad) * ad.T2.L.Z0 + Zprob.ad * ad.T2.L.Z1
    ad.T2.U <- (1 - Zprob.ad) * ad.T2.U.Z0 + Zprob.ad * ad.T2.U.Z1
    nd.T2.L <- (1 - Zprob.nd) * nd.T2.L.Z0 + Zprob.nd * nd.T2.L.Z1
    nd.T2.U <- (1 - Zprob.nd) * nd.T2.U.Z0 + Zprob.nd * nd.T2.U.Z1
    ad.T1.L <- (1 - Zprob.ad) * ad.T1.L.Z0 + Zprob.ad * ad.T1.L.Z1
    ad.T1.U <- (1 - Zprob.ad) * ad.T1.U.Z0 + Zprob.ad * ad.T1.U.Z1
    adjusted <- T
  }

    ret.list <- list(ad.T2.L = ad.T2.L, ad.T2.U = ad.T2.U, nd.T2.L = nd.T2.L, nd.T2.U = nd.T2.U,
                   ad.T1.L = ad.T1.L, ad.T1.U = ad.T1.U, adjusted = adjusted)
    ## bound for ad CDF
}
