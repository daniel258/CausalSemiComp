########################################################################
## Function for CausalSemiComp
### Estimation of various componenets for bounds and senstivity analysis

### The function returns the following:
########################################################################
CalcBounds <- function(Res, bounds = T, sens = F, sens.param = NULL, probs = T, RMST = T)
{
  all.times <- Res$all.times
  S2A0.all.times <- Res$S2A0.all.times
  S2A1.all.times <- Res$S2A0.all.times
  etaA0 <- Res$etaA0
  etaA1 <- Res$etaA1
  etasA0T2.le.t <- Res$etasA0T2.le.t
  etasA1T2.le.t <- Res$etasA1T2.le.t
  S1A0.all.times <- Res$S1A0.all.times
  S1A1.all.times <- Res$S1A1.all.times
  S1A0T1lT2.all.times <- res$S1A0T1lT2.all.times
  S1A1T1lT2.all.times <- res$S1A1T1lT2.all.times
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
  ret.list <- list(ad.T2.L = ad.T2.L, ad.T2.U = ad.T2.U, nd.T2.L = nd.T2.L, nd.T2.U = nd.T2.U,
                   ad.T1.L = ad.T1.L, ad.T1.U = ad.T1.U)
    ## bound for ad CDF
}
