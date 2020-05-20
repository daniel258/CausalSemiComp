BootCausalSC <- function(L = 0, T1, T2, delta1, delta2, A, all.times, B = 500)
{
  n.sample <- length(T1)

  for (b in 1:B)
  {
    indices <- sample(1:n.sample, size = n.sample, replace = T)
    T1.boot <- T1[indices]
    T2.boot <- T2[indices]
    delta1.boot <- delta1[indices]
    delta2.boot <- delta2[indices]
    A.boot <- A[indices]
    boot.etaA0 <- boot.etaA1 <- vector(length = n.sample)
    boot.S2A0 <- boot.S2A1 <- boot.S1A0 <- boot.S1A1 <-  matrix(nr = B, nc = n.times)
    boot.etasA0T2.le.t <- boot.etasA1T2.le.t <- boot.S1A0T1lT2.all <- boot.S1A1T1lT2.all <-
      matrix(nr = B, nc = n.times)
    if (L !=0) {
      L.boot <- L[indices]
      res.boot <- CausalSC(L = L.boot, T1 = T1.boot, T2 = T2.boot,
               delta1 = delta1.boot, delta2 = delta2.boot, A = A.boot,
               all.times = all.times)
    } else {
      res.boot <- CausalSC(T1 = T1.boot, T2 = T2.boot, delta1 = delta1.boot, delta2 = delta2.boot,
             A = A.boot, all.times = all.times)
    }
    boot.etaA0[b] <- res.boot$etaA0
    boot.etaA1[b] <- res.boot$etaA1
    boot.S2A0[b, ] <- res.boot$S2A0.all.times
    boot.S2A1[b, ] <- res.boot$S2A1.all.times
    boot.S1A0[b, ] <- res.boot$S1A0.all.times
    boot.S1A1[b, ] <- res.boot$S1A1.all.times
    boot.etasA0T2.le.t[b, ] <- res.boot$etasA0T2.le.t
    boot.etasA1T2.le.t[b, ] <- res.boot$etasA1T2.le.t
    boot.S1A0T1lT2.all[b, ] <- res.boot$S1A0T1lT2.all.times
    boot.S1A1T1lT2.all[b, ] <- res.boot$S1A1T1lT2.all.times
  }
  SE.etaA0 <- sd(boot.etaA0)
  SE.etaA1 <- sd(boot.etaA1)
  SE.S2A0 <- apply(boot.S2A0, 2, sd, na.rm = T)
  SE.S2A1 <- apply(boot.S2A1, 2, sd, na.rm = T)
  SE.S1A0 <- apply(boot.S1A0, 2, sd, na.rm = T)
  SE.S1A1 <- apply(boot.S1A1, 2, sd, na.rm = T)
  SE.etasA0T2.le.t <- apply(boot.etasA0T2.le.t, 2, sd, na.rm = T)
  SE.etasA1T2.le.t <- apply(boot.etasA1T2.le.t, 2, sd, na.rm = T)
  SE.S1A0T1lT2.all <- apply(boot.S1A0T1lT2.all, 2, sd, na.rm = T)
  SE.S1A1T1lT2.all <- apply(boot.S1A1T1lT2.all, 2, sd, na.rm = T)
  CI.np.etaA0 <- quantile(boot.etaA0, probs = c(0.025, 0.975), na.rm = T)
  CI.np.etaA1 <- quantile(boot.etaA1, probs = c(0.025, 0.975), na.rm = T)
  CI.np.S2A0 <- apply(boot.S2A0, 2, quantile, probs = c(0.025, 0.975), na.rm = T)
  CI.np.S2A1 <- apply(boot.S2A1, 2, quantile, probs = c(0.025, 0.975), na.rm = T)
  CI.np.S1A0 <- apply(boot.S1A0, 2, quantile, probs = c(0.025, 0.975), na.rm = T)
  CI.np.S1A1 <- apply(boot.S1A1, 2, quantile, probs = c(0.025, 0.975), na.rm = T)
  CI.np.etasA0T2.le.t <- apply(boot.etasA0T2.le.t, 2, quantile, probs = c(0.025, 0.975), na.rm = T)
  CI.np.etasA1T2.le.t <- apply(boot.etasA1T2.le.t, 2, quantile, probs = c(0.025, 0.975), na.rm = T)
  CI.np.S1A0T1lT2.all <- apply(boot.S1A0T1lT2.all, 2, quantile, probs = c(0.025, 0.975), na.rm = T)
  CI.np.S1A1T1lT2.all <- apply(boot.S1A1T1lT2.all, 2, quantile, probs = c(0.025, 0.975), na.rm = T)
  list.ret <- list(SE.etaA0 = SE.etaA0, SE.etaA1 = SE.etaA1, SE.S2A0 = SE.S2A0, SE.S2A1 = SE.S2A1,
                   SE.S1A0 = SE.S1A0, SE.S1A1 = SE.S1A1, SE.etasA0T2.le.t = SE.etasA0T2.le.t,
                   SE.etasA1T2.le.t = SE.etasA1T2.le.t, SE.S1A0T1lT2.all = SE.S1A0T1lT2.all,
                   SE.S1A1T1lT2.all = SE.S1A1T1lT2.all,
                   CI.np.etaA0 = CI.np.etaA0, CI.np.etaA1 = CI.np.etaA1, CI.np.S2A0 = CI.np.S2A0,
                   CI.np.S2A1 = CI.np.S2A1, CI.np.S1A0 = CI.np.S1A0, CI.np.S1A1 = CI.np.S1A1,
                   CI.np.etasA0T2.le.t = CI.np.etasA0T2.le.t, CI.np.etasA1T2.le.t = CI.np.etasA1T2.le.t,
                   CI.np.S1A0T1lT2.all = CI.np.S1A0T1lT2.all, CI.np.S1A1T1lT2.all = CI.np.S1A1T1lT2.all)
  return(list.ret)
}
