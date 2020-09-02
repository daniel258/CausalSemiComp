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

CausalSC <- function(L = 0, T1, T2, delta1, delta2, A, all.times)
{
  #if (!bounds & !sens) stop("Must specfiy either bounds=T or sens=T (or both)")
  #if (sens==T & is.null(thetas)) stop("Must specfiy thetas when sens==T")
  n.sample <- length(T1)
  n.times <- length(all.times)
  if  (!(length(L) %in% c(1, n.sample))) {
    stop("L should be either scalar or at same length as T1")
  }
  # Estimate S_{2|A=a}(t)
  if(length(L)==1)
  {
    S2A0 <- survival::survfit(survival::Surv(T2, delta2) ~ 1, subset = A==0)
    S2A1 <- survival::survfit(survival::Surv(T2, delta2) ~ 1, subset = A==1)
  } else {
    S2A0 <- survival::survfit(survival::Surv(L, T2, delta2) ~ 1, subset = A==0)
    S2A1 <- survival::survfit(survival::Surv(L, T2, delta2) ~ 1, subset = A==1)
  }
  S2A0.surv <- S2A0$surv
  S2A1.surv <- S2A1$surv
  dS2A0 <-  diff(c(1, S2A0.surv))
  dS2A1 <-  diff(c(1, S2A1.surv))
  dS2A0.times <-  S2A0$time
  dS2A1.times <-  S2A1$time
  n.times.A0 <- length(dS2A0.times)
  n.times.A1 <- length(dS2A1.times)
  # Calculate Beran's smoothed estimator for S_{1|A=a, T_2=t_t}
  T1A0dead <- T1[A==0 & delta2==1]
  T1A1dead <- T1[A==1 & delta2==1]
  delta1A0dead <- delta1[A==0 & delta2==1]
  delta1A1dead <- delta1[A==1 & delta2==1]
  T2A0dead <- T2[A==0 & delta2==1]
  T2A1dead <- T2[A==1 & delta2==1]
  if(length(L)==1)
  {
    fit.eta.A0 <- prodlim::prodlim(survival::Surv(T1A0dead, delta1A0dead) ~ T2A0dead)
    fit.eta.A1 <- prodlim::prodlim(survival::Surv(T1A1dead, delta1A1dead) ~ T2A1dead)
  } else {
    LA0dead <- L[A==0 & delta2==1]
    LA1dead <- L[A==1 & delta2==1]
    fit.eta.A0 <- prodlim::prodlim(survival::Surv(LA0dead, T1A0dead, delta1A0dead) ~ T2A0dead)
    fit.eta.A1 <- prodlim::prodlim(survival::Surv(LA1dead, T1A1dead, delta1A1dead) ~ T2A1dead)
  }
  # Calculate S_{1|A=0, T_2=t_t} only for t for which dS2A0 > 0 (and the same for A=1 with dS2A1 > 0)
  newdata.A0 <- data.frame(T2A0dead = S2A0$time)
  newdata.A1 <- data.frame(T2A1dead = S2A1$time)
  sorted.T2A0dead.times <- sort(newdata.A0$T2A0dead)
  sorted.T2A1dead.times <- sort(newdata.A1$T2A1dead)
  Beran.curves.A0.T1 <- predict(fit.eta.A0, newdata = newdata.A0, times = sorted.T2A0dead.times)
  Beran.curves.A1.T1 <- predict(fit.eta.A1, newdata = newdata.A1, times = sorted.T2A1dead.times)
  Beran.curves.A0.T1.mat <- do.call(rbind, Beran.curves.A0.T1)
  Beran.curves.A1.T1.mat <- do.call(rbind, Beran.curves.A1.T1)
  # Change NA after no more events into the last estimated value
  Beran.curves.A0.T1.mat <- t(apply(Beran.curves.A0.T1.mat, 1, RepNAmin))
  Beran.curves.A1.T1.mat <- t(apply(Beran.curves.A1.T1.mat, 1, RepNAmin))
  diag.beran.mat.A0.T1 <- diag(Beran.curves.A0.T1.mat)
  diag.beran.mat.A1.T1 <- diag(Beran.curves.A1.T1.mat)
  #### Calculate eta_{A=a} ####
  S2A0.tau <- S2A0$surv[length(S2A0$surv)]
  S2A1.tau <- S2A1$surv[length(S2A1$surv)]
  etaA0 <- -sum((1 - diag.beran.mat.A0.T1) * dS2A0)
  etaA1 <- -sum((1 - diag.beran.mat.A1.T1) * dS2A1)
  #etaA0 <- -sum((1 - diag.beran.mat.A0.T1) * dS2A0)/(1 - S2A0.tau)
  #etaA1 <- -sum((1 - diag.beran.mat.A1.T1) * dS2A1)/(1 - S2A1.tau)
  ## Define a function that returns S_{1|A=a,T_2=s}(t)
   S1AtT2s <- function(my.t, my.s, A) {
     if (my.s < my.t)
     {
       my.t <- my.s
     }
     if (A==0) {
     s.place <- findInterval(my.s, newdata.A0$T2A0dead)
     t.place <- findInterval(my.t, newdata.A0$T2A0dead)
     if (t.place==0 | s.place==0) {
       S.res <- 1
     } else {
       t1.curve.s <- Beran.curves.A0.T1.mat[s.place, ]
       if (t.place==length(newdata.A0$T2A0dead)) {
         S.res <-  t1.curve.s[length(newdata.A0$T2A0dead)]
     } else {
       S.res <- t1.curve.s[t.place]
     }}}
     if (A==1) {
       s.place <- findInterval(my.s, newdata.A1$T2A1dead)
       t.place <- findInterval(my.t, newdata.A1$T2A1dead)
       if (t.place==0 | s.place==0) {
         S.res <- 1
       } else {
         t1.curve.s <- Beran.curves.A1.T1.mat[s.place,]
         if (t.place==length(newdata.A1$T2A0dead)) {
           S.res <-  t1.curve.s[length(newdata.A1$T2A1dead)]
         } else {
           S.res <- t1.curve.s[t.place]
         }}}
   return(S.res)
   }
   #### Calculate S_{1|A=a}(t) ####
   S1A0.all.times <- S1A1.all.times <- vector(length = n.times)
   S1sA0t <- matrix(nr = n.times, nc = n.times.A0)
   S1sA1t <- matrix(nr = n.times, nc = n.times.A1)
   for (j in 1:n.times)
    {
    t.now <- all.times[j]
    A0.t.now <- findInterval(t.now, dS2A0.times)
    A1.t.now <- findInterval(t.now, dS2A1.times)
    #t.now.before <- all(dS2A0.times > t.now)
     if(A0.t.now==0) {
       S1A0.all.times[j] <- 1
     } else {
       S1sA0t[j, ] <-  sapply(dS2A0.times, S1AtT2s, my.t = t.now, A = 0)
       S1A0.all.times[j] <- -sum(S1sA0t[j, ] * dS2A0)
       #S1A0.all.times[j] <- -sum(S1sA0t[j, ] * dS2A0) / (1 - S2A0.tau)
     }
     if (A1.t.now==0) {
       S1A1.all.times[j] <- 1
     } else {
       S1sA1t[j, ] <-  sapply(dS2A1.times, S1AtT2s, my.t = t.now, A = 1)
       S1A1.all.times[j] <- -sum(S1sA1t[j, ] * dS2A1)
     #  S1A1.all.times[j] <- -sum(S1sA1t[j, ] * dS2A1) / (1 - S2A1.tau)
     }}
    #### eta_{A=a,T_2\le t} and S_{1|A=a,T1 \le T_2} ####
    etasA0T2.le.t <- etasA1T2.le.t <- vector(length = n.times)
    S1A0T1lT2.all.times <- S1A1T1lT2.all.times <- vector(length = n.times)
    #### Calculate eta_{A=a,T_2\le t} ####
    for (j in 1:n.times)
    {
      t.now <- all.times[j]
      A0.t.now <- findInterval(t.now, dS2A0.times)
      A1.t.now <- findInterval(t.now, dS2A1.times)
    #A0.t.now <- findInterval(t.now, dS2A0.times)
    if (A0.t.now==0) {
      S2A0.now <- 1
      } else if (A0.t.now==n.times.A0) {
        S2A0.now <- S2A0.surv[n.times.A0]
        } else {
      S2A0.now <- S2A0.surv[A0.t.now]
        }
    #A1.t.now <- findInterval(t.now, dS2A1.times)
    if (A1.t.now==0) {
      S2A1.now <- 1
    } else if (A1.t.now==n.times.A1) {
      S2A1.now <- S2A1.surv[n.times.A1]
    } else {
      S2A1.now <- S2A1.surv[A1.t.now]
    }
    etasA0T2.le.t[j] <- -sum((1 - diag.beran.mat.A0.T1[dS2A0.times <= t.now]) * dS2A0[dS2A0.times <= t.now]) / (1 - S2A0.now)
    etasA1T2.le.t[j] <- -sum((1 - diag.beran.mat.A1.T1[dS2A1.times <= t.now]) * dS2A1[dS2A1.times <= t.now]) / (1 - S2A1.now)
    if (A0.t.now==0) {
      S1A0T1lT2.all.times[j] <- 1
    } else {
    S1A0T1lT2.all.times[j] <-  -sum((Beran.curves.A0.T1.mat[dS2A0.times > t.now, A0.t.now] -
                                       diag.beran.mat.A0.T1[dS2A0.times > t.now]) * dS2A0[dS2A0.times > t.now]) #/ (1 - S2A0.tau)
    }
    if (A1.t.now==0) {
      S1A1T1lT2.all.times[j] <- 1
    } else {
      S1A1T1lT2.all.times[j] <-  -sum((Beran.curves.A1.T1.mat[dS2A1.times > t.now, A1.t.now] -
                                         diag.beran.mat.A1.T1[dS2A1.times > t.now]) * dS2A1[dS2A1.times > t.now]) #/ (1 - S2A1.tau)
    }}
    S1A0T1lT2.all.times <- S1A0T1lT2.all.times/etaA0
    S1A1T1lT2.all.times <- S1A1T1lT2.all.times/etaA1

    S2A0Func <- stepfun(x = S2A0$time , y = c(1, S2A0$surv), right = F)
    S2A1Func <- stepfun(x = S2A1$time , y = c(1, S2A1$surv), right = F)
    S2A0.all.times <- sapply(all.times, S2A0Func)
    S2A1.all.times <- sapply(all.times, S2A1Func)
    surv.comp <- list(S2A0 = S2A0, S2A1 = S2A1,
                      S2A0.all.times = S2A0.all.times,  S2A1.all.times = S2A1.all.times,
                      etaA0 = etaA0, etaA1 = etaA1,
                      etasA0T2.le.t = etasA0T2.le.t, etasA1T2.le.t = etasA1T2.le.t,
                      S1A0.all.times = S1A0.all.times, S1A1.all.times = S1A1.all.times,
                      S1A0T1lT2.all.times = S1A0T1lT2.all.times, S1A1T1lT2.all.times = S1A1T1lT2.all.times,
                      all.times = all.times)
  return(surv.comp)
}

