WrapCausalSCboot <- function(data, i = i, all.times, Ltime, out.times)
{
  if (Ltime==0)
  {
  res <- CausalSC(T1 = data$T1[i], T2 = data$T2[i],
                  delta1 = data$delta1[i], delta2 = data$delta2[i], A = data$A[i],
                  all.times = all.times) }
  else {
    res <- CausalSC(T1 = data$T1[i], T2 = data$T2[i],
                    delta1 = data$delta1[i], delta2 = data$delta2[i], A = data$A[i], L = data$L[i],
                    all.times = all.times)
  }
  ind.times <- all.times %in% out.times
  res.out <- c(res$etaA0, res$etaA1, res$S2A0.all.times[ind.times], res$S2A1.all.times[ind.times],
               res$S1A0.all.times[ind.times], res$S1A1.all.times[ind.times],
               res$etasA0T2.le.t[ind.times], res$etasA1T2.le.t[ind.times],
               res$S1A0T1lT2.all.times[ind.times], res$S1A1T1lT2.all.times[ind.times])
  return(res.out)
}
