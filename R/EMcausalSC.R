#' @title Estimation of components of the bounds for population stratification effects and stratum proportion
#' @description The function implements the proposed methodology in Nevo and Gorfine (2021+, Biostatistics)
#' @details bla bla bla
#' @param rho
#' @param tau
#' @param n.gamma.vals
#' @param n.sample.pers
#' @param population
#' @param Xnames
#' @param data
#' @param res
#' @return S_{2|A=a}(t), for a=0,1  (names: S2A0 and S2A1)
#  \eta_{A=a}, for a=0,1  (names: etaA0 and etaA1)
#  \eta_{A=a, T_2 \le t} for a=0,1 (names: etasA0T2.le.t and etasA1T2.le.t)
#  S_{1|A=a}(t) for a=0,1 (names: S1A0.all.times and S1A1.all.times)
#  all.times
#' @examples
#' \dontrun{
#' set.seed(314)
#' times <- seq(1, 15, 1)
#' alpha.nt <- LongitSemiComp:::logit(dchisq(times,3, ncp =5)/2 + 0.025)
#' alpha.t <- LongitSemiComp:::logit(times*(0.075/10)  - 0.0005*(times/20)^2  + 0.05)
#' alpha.or <- 0.15 - times/10 + 0.75*(times/10)^2 + 0.3*(times/20)^3
#' plot(x = times, y= exp(alpha.or))
#' plot(x = times, y= LongitSemiComp:::expit(alpha.nt))
#' plot(x = times, y= LongitSemiComp:::expit(alpha.t))
#' beta.nt <- log(c(0.7, 3))
#' beta.t <- log(c(0.5, 1))
#' beta.or <- log(c(1.4, 1))
#' beta.y <- log(1.4)
#' my.data <- SimLongitData(n.sample = 2000, times = times,  beta.y,
#'                          alpha.nt, alpha.t, alpha.or,
#'                          beta.nt, beta.t, beta.or)
#' longit.data <- my.data[-1]
#' X <- as.data.frame(my.data[1])
#' LongitSC(longit.data = longit.data, times = times,  formula.NT =  ~ X.1 + X.2,
#'          formula.T =  ~ X.1 + X.2,
#'          formula.OR = ~ X.1 + X.2,
#'          data = X, epsOR = 10^(-10),
#'          knots = 5, lambda = 1)
#' }
#' @author Daniel Nevo
#' @importFrom stats as.formula optim optimize
#' @export
EMcausalSC <- function(data, Xnames, Lname = NULL,
                       max.iter = 10000, w = NULL, eps.conv = 0.0001, init.thetas = c(1, 1))
{
  n.sample <- nrow(data)
  delta1 <- data$delta1
  delta2 <- data$delta2
  A <- data$A
  #### If weights were supplied, normalize by A value ###
  if (!is.null(w)) {
    if (length(w)!=n.sample) {stop("w should be of the same length as the sample size")}
    w[A==0] <- w[A==0]/mean(w[A==0])
    w[A==1] <- w[A==1]/mean(w[A==1])
  data$w <- w  # add the weights to the data
  data <- data[w > 0, ]
  A <- A[w >0]
  delta1 <- delta1[w > 0]
  delta2 <- delta2[w > 0]
  n.sample <- length(A)
  }
  delta1A0 <- delta1[A==0]
  delta1A1 <- delta1[A==1]
  delta2A0 <- delta2[A==0]
  delta2A1 <- delta2[A==1]
  data$delta2only <- data$delta2*(1 - data$delta1)
  # data includes T1, T2, delta1, delta2, A, X
  ## This is used for calculating
  data$log.gamma = 0
  ### Create data subsets for the different Cox models
  new.data.A0 <- data %>% dplyr::filter(A==0)
  new.data.A1 <- data %>% dplyr::filter(A==1)
  new.data.A0T12 <- data %>% dplyr::filter(A==0 & delta1==1)
  new.data.A1T12 <- data %>% dplyr::filter(A==1 & delta1==1)
  XmatA0 <- as.matrix(dplyr::select(new.data.A0, Xnames))
  XmatA1 <- as.matrix(dplyr::select(new.data.A1, Xnames))
  XmatA0T12 <- as.matrix(dplyr::select(new.data.A0T12, Xnames))
  XmatA1T12 <- as.matrix(dplyr::select(new.data.A1T12, Xnames))
  m.X.A0 <- apply(XmatA0, 2, mean)
  m.X.A1 <- apply(XmatA1, 2, mean)
  m.X.A0T12 <- apply(XmatA0T12, 2, mean)
  m.X.A1T12 <- apply(XmatA1T12, 2, mean)
  ## Add psuedo observations with no truncation and mean X, logGamma for survival function calculations:
  data.predict.A0 <- new.data.A0[1, ]
  data.predict.A1 <- new.data.A1[1, ]
  data.predict.A0T12 <- new.data.A0T12[1, ]
  data.predict.A1T12 <- new.data.A1T12[1, ]
  # No truncation
  data.predict.A0[, names(data.predict.A0)==Lname] <- data.predict.A1[, names(data.predict.A1)==Lname] <- 0
  data.predict.A0T12[, names(data.predict.A0T12)==Lname] <- 0
  data.predict.A1T12[, names(data.predict.A1T12)==Lname] <- 0
  # Mean X values
  data.predict.A0[, names(data.predict.A0) %in% Xnames]  <- 0
  data.predict.A1[, names(data.predict.A1) %in% Xnames]  <- 0
  data.predict.A0T12[, names(data.predict.A0T12) %in% Xnames]  <- 0
  data.predict.A1T12[, names(data.predict.A1T12) %in% Xnames]  <- 0
  #### Formulas for the analysis
  if(is.null(Lname))
  {
  formula01 <- as.formula(paste0("Surv(T1, delta1) ~", paste0(Xnames, collapse = "+"),
                                 " + offset(log.gamma)"))
  formula02 <- as.formula(paste0("Surv(T1, delta2only) ~", paste0(Xnames, collapse = "+"),
                                 " + offset(log.gamma)"))
  formula12 <- as.formula(paste0("Surv(T1, T2, delta2) ~",
                                 paste0(Xnames, collapse= "+") ,
                                 " + offset(log.gamma)"))
  } else {
  formula01 <- as.formula(paste0("Surv(", Lname, ", T1, delta1) ~", paste0(Xnames, collapse = "+"),
                                 " + offset(log.gamma)"))
  formula02 <- as.formula(paste0("Surv(", Lname, ",T1, delta2only) ~", paste0(Xnames, collapse = "+"),
                                   " + offset(log.gamma)"))
  formula12 <- as.formula(paste0("Surv(T1, T2, delta2) ~",
                                   paste0(Xnames, collapse= "+") ,
                                   " + offset(log.gamma)"))
  }


  ########  Initial parameter values ########
  if (is.null(w))
  {
    fit.a0.01 <- survival::coxph(formula01, data = new.data.A0)
    fit.a0.02 <- survival::coxph(formula02, data = new.data.A0)
    fit.a0.12 <- survival::coxph(formula12, data = new.data.A0T12)
    fit.a1.01 <- survival::coxph(formula01, data = new.data.A1)
    fit.a1.02 <- survival::coxph(formula02, data = new.data.A1)
    fit.a1.12 <- survival::coxph(formula12, data = new.data.A1T12)
  } else {
    fit.a0.01 <- survival::coxph(formula01, data = new.data.A0, weights = w)
    fit.a0.02 <- survival::coxph(formula02, data = new.data.A0, weights = w)
    fit.a0.12 <- survival::coxph(formula12, data = new.data.A0T12, weights = w)
    fit.a1.01 <- survival::coxph(formula01, data = new.data.A1, weights = w)
    fit.a1.02 <- survival::coxph(formula02, data = new.data.A1, weights = w)
    fit.a1.12 <- survival::coxph(formula12, data = new.data.A1T12, weights = w)
  }
  est.beta.a0.01 <- coef(fit.a0.01)
  est.beta.a1.01 <- coef(fit.a1.01)
  est.beta.a0.02 <- coef(fit.a0.02)
  est.beta.a1.02 <- coef(fit.a1.02)
  est.beta.a0.12  <- coef(fit.a0.12)
  est.beta.a1.12  <- coef(fit.a1.12)
  old.betas <- new.betas <- naive.betas <-  c(est.beta.a0.01, est.beta.a0.02, est.beta.a0.12,
                                              est.beta.a1.01, est.beta.a1.02, est.beta.a1.12)
  old.thetas <- new.thetas <- init.thetas
  # data.predict.A0$log.gamma <- mean(new.data.A0$log.gamma)
  # data.predict.A1$log.gamma <- mean(new.data.A1$log.gamma)
  # data.predict.A0T12$log.gamma <- mean(new.data.A0T12$log.gamma)
  # data.predict.A1T12$log.gamma <- mean(new.data.A1T12$log.gamma)
  s.fit.a0.1 <- survival::survfit(fit.a0.01, censor = FALSE, se.fit = FALSE,  newdata = data.predict.A0)
  s.fit.a0.2 <- survival::survfit(fit.a0.02, censor = FALSE, se.fit = FALSE,  newdata = data.predict.A0)
  s.fit.a0.12 <- survival::survfit(fit.a0.12, censor = FALSE, se.fit = FALSE, newdata = data.predict.A0T12)
  s.fit.a1.1 <- survival::survfit(fit.a1.01, censor = FALSE, se.fit = FALSE,  newdata = data.predict.A1)
  s.fit.a1.2 <- survival::survfit(fit.a1.02,  censor = FALSE, se.fit = FALSE, newdata = data.predict.A1)
  s.fit.a1.12 <- survival::survfit(fit.a1.12, censor = FALSE, se.fit = FALSE, newdata = data.predict.A1T12)
  step.A0T1 <- stepfun(x = s.fit.a0.1$time,   y = c(0, s.fit.a0.1$cumhaz))
  step.A0T2 <- stepfun(x = s.fit.a0.2$time,   y = c(0, s.fit.a0.2$cumhaz))
  step.A0T12 <- stepfun(x = s.fit.a0.12$time, y = c(0, s.fit.a0.12$cumhaz))
  step.A1T1 <- stepfun(x = s.fit.a1.1$time,   y = c(0, s.fit.a1.1$cumhaz))
  step.A1T2 <- stepfun(x = s.fit.a1.2$time,   y = c(0, s.fit.a1.2$cumhaz))
  step.A1T12 <- stepfun(x = s.fit.a1.12$time, y = c(0, s.fit.a1.12$cumhaz))
  # step.A0T1 <- stepfun(x = s.fit.a0.1$time,
  #                      y = c(0, s.fit.a0.1$cumhaz*exp(-sum(est.beta.a0.01*m.X.A0)
  #                                                     -mean(new.data.A0$log.gamma))))
  # step.A0T2 <- stepfun(x = s.fit.a0.2$time,
  #                      y = c(0, s.fit.a0.2$cumhaz*exp(-sum(est.beta.a0.02*m.X.A0)
  #                                                     - mean(new.data.A0$log.gamma))))
  # step.A0T12 <- stepfun(x = s.fit.a0.12$time,
  #                       y = c(0, s.fit.a0.12$cumhaz*exp(-sum(est.beta.a0.12*m.X.A0T12)
  #                                                       - mean(new.data.A0T12$log.gamma))))
  # step.A1T1 <- stepfun(x = s.fit.a1.1$time,
  #                      y = c(0, s.fit.a1.1$cumhaz*exp(-sum(est.beta.a1.01*m.X.A1)
  #                                                     -mean(new.data.A1$log.gamma))))
  # step.A1T2 <- stepfun(x = s.fit.a1.2$time,
  #                      y = c(0, s.fit.a1.2$cumhaz*exp(-sum(est.beta.a1.02*m.X.A1)
  #                                                     - mean(new.data.A1$log.gamma))))
  # step.A1T12 <- stepfun(x = s.fit.a1.12$time,
  #                       y = c(0, s.fit.a1.12$cumhaz*exp(-sum(est.beta.a1.12*m.X.A1T12)
  #                                                       - mean(new.data.A1T12$log.gamma))))
  ################################################################################
  ################################################################################################
  # ##### Calculate per-person cumulative hazards (with gamma=0)
  # ################################################################################################

  H.a0.01 <- step.A0T1(new.data.A0$T1)        * exp(XmatA0%*%est.beta.a0.01)
  H.a0.02 <- step.A0T2(new.data.A0$T1)        * exp(XmatA0%*%est.beta.a0.02)
  H.a0.12 <- step.A0T12(new.data.A0T12$T2)    * exp(XmatA0T12%*%est.beta.a0.12)
  H.a0.12.T1 <- step.A0T12(new.data.A0T12$T1) * exp(XmatA0T12%*%est.beta.a0.12)
  H.a1.01 <- step.A1T1(new.data.A1$T1)        * exp(XmatA1%*%est.beta.a1.01)
  H.a1.02 <- step.A1T2(new.data.A1$T1)        * exp(XmatA1%*%est.beta.a1.02)
  H.a1.12 <- step.A1T12(new.data.A1T12$T2)    * exp(XmatA1T12%*%est.beta.a1.12)
  H.a1.12.T1 <- step.A1T12(new.data.A1T12$T1) * exp(XmatA1T12%*%est.beta.a1.12)
  if (!is.null(Lname))
  {
    L.a0.01 <- step.A0T1(new.data.A0[, names(new.data.A0)==Lname]) * exp(XmatA0%*%est.beta.a0.01)
    L.a0.02 <- step.A0T2(new.data.A0[, names(new.data.A0)==Lname]) * exp(XmatA0%*%est.beta.a0.02)
    L.a1.01 <- step.A1T1(new.data.A1[, names(new.data.A1)==Lname]) * exp(XmatA1%*%est.beta.a1.01)
    L.a1.02 <- step.A1T2(new.data.A1[, names(new.data.A1)==Lname]) * exp(XmatA1%*%est.beta.a1.02)
  }
  iter <- cond <- 0

  E.gamma <- E.log.gamma <- s.i <- vector(length = n.sample)
  ### Finally, the EM loop
  while(cond==0 & iter < max.iter)
  {
    iter <- iter + 1
    # Daniel::CatIndex(iter)
    # Daniel::CatIndex(round(new.thetas, 2))
    ##### E-step
    #####  Per-person posterior distriubtion parametrs
    if (is.null(Lname))
    {
      s.i[A==0] <- H.a0.01 + H.a0.02
      s.i[A==1] <- H.a1.01 + H.a1.02
    }
    else {
      s.i[A==0] <- H.a0.01 + H.a0.02 - L.a0.01 - L.a0.02
      s.i[A==1] <- H.a1.01 + H.a1.02 - L.a1.01 - L.a1.02
    }
    s.i[A==0 & delta1==1] <- s.i[A==0 & delta1==1] + H.a0.12 - H.a0.12.T1
    s.i[A==1 & delta1==1] <- s.i[A==1 & delta1==1] + H.a1.12 - H.a1.12.T1
    E.gamma[A==0] <- (1/new.thetas[1] + delta1A0 + delta2A0) / (1/new.thetas[1] + s.i[A==0])
    E.gamma[A==1] <- (1/new.thetas[2] + delta1A1 + delta2A1) / (1/new.thetas[2] + s.i[A==1])
    E.log.gamma[A==0] <- digamma(1/new.thetas[1] + delta1A0 + delta2A0) - log(1/new.thetas[1] + s.i[A==0])
    E.log.gamma[A==1] <- digamma(1/new.thetas[2] + delta1A1 + delta2A1) - log(1/new.thetas[2] + s.i[A==1])
    log.E.gamma <- log(E.gamma)
    new.data.A0$log.gamma <- log.E.gamma[A==0]
    new.data.A1$log.gamma <- log.E.gamma[A==1]
    new.data.A0T12$log.gamma <- log.E.gamma[A==0 & delta1==1]
    new.data.A1T12$log.gamma <- log.E.gamma[A==1 & delta1==1]
    ###############################################################################################
    ##### M-step
    #### Conditonially on gamma, fit illness-death PH models
    ###############################################################################################
    if (is.null(w))
    {
      fit.a0.01 <- survival::coxph(formula01, data = new.data.A0)
      fit.a0.02 <- survival::coxph(formula02, data = new.data.A0)
      fit.a0.12 <- survival::coxph(formula12, data = new.data.A0T12)
      fit.a1.01 <- survival::coxph(formula01, data = new.data.A1)
      fit.a1.02 <- survival::coxph(formula02, data = new.data.A1)
      fit.a1.12 <- survival::coxph(formula12, data = new.data.A1T12)
    } else {
      fit.a0.01 <- survival::coxph(formula01, data = new.data.A0, weights = w)
      fit.a0.02 <- survival::coxph(formula02, data = new.data.A0, weights = w)
      fit.a0.12 <- survival::coxph(formula12, data = new.data.A0T12, weights = w)
      fit.a1.01 <- survival::coxph(formula01, data = new.data.A1, weights = w)
      fit.a1.02 <- survival::coxph(formula02, data = new.data.A1, weights = w)
      fit.a1.12 <- survival::coxph(formula12, data = new.data.A1T12, weights = w)
    }
    est.beta.a0.01 <- coef(fit.a0.01)
    est.beta.a1.01 <- coef(fit.a1.01)
    est.beta.a0.02 <- coef(fit.a0.02)
    est.beta.a1.02 <- coef(fit.a1.02)
    est.beta.a0.12  <- coef(fit.a0.12)
    est.beta.a1.12  <- coef(fit.a1.12)
    new.betas <- c(est.beta.a0.01, est.beta.a0.02, est.beta.a0.12,
                   est.beta.a1.01, est.beta.a1.02, est.beta.a1.12)
    ################################################################################################
    ##### Create step functions from all baseline hazard estimators ########
    ################################################################################################
    # data.predict.A0$log.gamma <- mean(new.data.A0$log.gamma)
    # data.predict.A1$log.gamma <- mean(new.data.A1$log.gamma)
    # data.predict.A0T12$log.gamma <- mean(new.data.A0T12$log.gamma)
    # data.predict.A1T12$log.gamma <- mean(new.data.A1T12$log.gamma)
    s.fit.a0.1 <- survival::survfit(fit.a0.01, censor = FALSE, se.fit = FALSE,  newdata = data.predict.A0)
    s.fit.a0.2 <- survival::survfit(fit.a0.02, censor = FALSE, se.fit = FALSE,  newdata = data.predict.A0)
    s.fit.a0.12 <- survival::survfit(fit.a0.12, censor = FALSE, se.fit = FALSE, newdata = data.predict.A0T12)
    s.fit.a1.1 <- survival::survfit(fit.a1.01, censor = FALSE, se.fit = FALSE,  newdata = data.predict.A1)
    s.fit.a1.2 <- survival::survfit(fit.a1.02,  censor = FALSE, se.fit = FALSE, newdata = data.predict.A1)
    s.fit.a1.12 <- survival::survfit(fit.a1.12, censor = FALSE, se.fit = FALSE, newdata = data.predict.A1T12)
    step.A0T1 <- stepfun(x = s.fit.a0.1$time,   y = c(0, s.fit.a0.1$cumhaz))
    step.A0T2 <- stepfun(x = s.fit.a0.2$time,   y = c(0, s.fit.a0.2$cumhaz))
    step.A0T12 <- stepfun(x = s.fit.a0.12$time, y = c(0, s.fit.a0.12$cumhaz))
    step.A1T1 <- stepfun(x = s.fit.a1.1$time,   y = c(0, s.fit.a1.1$cumhaz))
    step.A1T2 <- stepfun(x = s.fit.a1.2$time,   y = c(0, s.fit.a1.2$cumhaz))
    step.A1T12 <- stepfun(x = s.fit.a1.12$time, y = c(0, s.fit.a1.12$cumhaz))
    # step.A0T1 <- stepfun(x = s.fit.a0.1$time,
    #                      y = c(0, s.fit.a0.1$cumhaz*exp(- sum(est.beta.a0.01*m.X.A0)
    #                                                     - mean(new.data.A0$log.gamma))))
    # step.A0T2 <- stepfun(x = s.fit.a0.2$time,
    #                      y = c(0, s.fit.a0.2$cumhaz*exp(- sum(est.beta.a0.02*m.X.A0)
    #                                                     - mean(new.data.A0$log.gamma))))
    #
    # step.A0T12 <- stepfun(x = s.fit.a0.12$time,
    #                       y = c(0, s.fit.a0.12$cumhaz*exp(- sum(est.beta.a0.12*m.X.A0T12)
    #                                                       - mean(new.data.A0T12$log.gamma))))
    # step.A1T1 <- stepfun(x = s.fit.a1.1$time,
    #                      y = c(0, s.fit.a1.1$cumhaz*exp(- sum(est.beta.a1.01*m.X.A1)
    #                                                     - mean(new.data.A1$log.gamma))))
    # step.A1T2 <- stepfun(x = s.fit.a1.2$time,
    #                      y = c(0, s.fit.a1.2$cumhaz*exp(- sum(est.beta.a1.02*m.X.A1)
    #                                                     - mean(new.data.A1$log.gamma))))
    # step.A1T12 <- stepfun(x = s.fit.a1.12$time,
    #                       y = c(0, s.fit.a1.12$cumhaz*exp(-sum(est.beta.a1.12*m.X.A1T12)
    #                                                       - mean(new.data.A1T12$log.gamma))))
    # Per-person cumulative hazards (with gamma=0) for "posterior distribution"
    ##################################################################################################
    H.a0.01 <- step.A0T1(new.data.A0$T1) * exp(XmatA0%*%est.beta.a0.01)
    H.a0.02 <- step.A0T2(new.data.A0$T1) * exp(XmatA0%*%est.beta.a0.02)
    H.a0.12 <- step.A0T12(new.data.A0T12$T2) * exp(XmatA0T12%*%est.beta.a0.12)
    H.a0.12.T1 <- step.A0T12(new.data.A0T12$T1) * exp(XmatA0T12%*%est.beta.a0.12)
    H.a1.01 <- step.A1T1(new.data.A1$T1) * exp(XmatA1%*%est.beta.a1.01)
    H.a1.02 <- step.A1T2(new.data.A1$T1) * exp(XmatA1%*%est.beta.a1.02)
    H.a1.12 <- step.A1T12(new.data.A1T12$T2) * exp(XmatA1T12%*%est.beta.a1.12)
    H.a1.12.T1 <- step.A1T12(new.data.A1T12$T1) * exp(XmatA1T12%*%est.beta.a1.12)

    if (!is.null(Lname))
    {
      L.a0.01 <- step.A0T1(new.data.A0[, names(new.data.A0)==Lname]) * exp(XmatA0%*%est.beta.a0.01)
      L.a0.02 <- step.A0T2(new.data.A0[, names(new.data.A0)==Lname]) * exp(XmatA0%*%est.beta.a0.02)
      L.a1.01 <- step.A1T1(new.data.A1[, names(new.data.A1)==Lname]) * exp(XmatA1%*%est.beta.a1.01)
      L.a1.02 <- step.A1T2(new.data.A1[, names(new.data.A1)==Lname]) * exp(XmatA1%*%est.beta.a1.02)
    }
    ##### New theta values #####
    if (is.null(w))
    {
      new.thetas[1] <- optimize(f = Eloglik, interval = c(0.0001, 30),
                                delta1 = delta1A0, delta2 = delta2A0,
                                E.gamma = E.gamma[A==0],
                                E.log.gamma = E.log.gamma[A==0],
                                maximum = T)$maximum
      new.thetas[2] <- optimize(f = Eloglik, interval = c(0.0001, 15),
                                delta1 = delta1A1, delta2 = delta2A1,
                                E.gamma = E.gamma[A==1],
                                E.log.gamma = E.log.gamma[A==1],
                                maximum = T)$maximum
    } else {
        new.thetas[1] <- optimize(f = Eloglik, interval = c(0.0001, 15),
                              delta1 = delta1A0, delta2 = delta2A0,
                              E.gamma = E.gamma[A==0],
                              E.log.gamma = E.log.gamma[A==0], w = w[A==0],
                              maximum = T)$maximum
        new.thetas[2] <- optimize(f = Eloglik, interval = c(0.0001, 15),
                              delta1 = delta1A1, delta2 = delta2A1,
                              E.gamma = E.gamma[A==1],
                              E.log.gamma = E.log.gamma[A==1], w = w[A==1],
                              maximum = T)$maximum
    }
    if (max(abs(c(new.thetas - old.thetas, new.betas - old.betas))) < eps.conv) {cond <- 1}
    old.betas <- new.betas
    old.thetas <- new.thetas
  }
  fit.list <- list(fit.a0.01 = fit.a0.01, fit.a0.02 = fit.a0.02, fit.a0.12 = fit.a0.12,
                   fit.a1.01 = fit.a1.01, fit.a1.02 = fit.a1.02, fit.a1.12 = fit.a1.12)
  H.step.funcs <- list(step.A0T1 = step.A0T1, step.A0T2 = step.A0T2,
                       step.A0T12 = step.A0T12,  step.A1T1 = step.A1T1,
                       step.A1T2 = step.A1T2, step.A1T12 = step.A1T12)
  list.out <- list(betas = new.betas, thetas = new.thetas, naive.betas = naive.betas,
                   fit.list = fit.list, H.step.funcs = H.step.funcs, iter = iter, E.gamma = E.gamma, mean.A = mean(A))
  return(list.out)
}
