####################################################################################
####################################################################################
# CausalSemiComp
# EM function for frailty SCM
####################################################################################
EMcausalSC <- function(data, Xnames, max.iter = 10000)
{
  n.sample <- nrow(data)
  Ltime <- "L" %in% colnames(data)
  delta1 <- data$delta1
  delta2 <- data$delta2
  A <- data$A
  delta1A0 <- delta1[A==0]
  delta1A1 <- delta1[A==1]
  delta2A0 <- delta2[A==0]
  delta2A1 <- delta2[A==1]
  data$delta2only <- data$delta2*(1 - data$delta1)
  # data includes T1, T2, delta1, delta2, A, X
  ## This is used for calculating
  data$log.gamma = 0
  ### Create data subsets for the different Cox models
  new.data.A0 <- data.predictA0 <- data %>% dplyr::filter(A==0)
  new.data.A1 <- data.predictA1 <- data %>% dplyr::filter(A==1)
  new.data.A0T12 <- data.predictA0T12 <- data %>% dplyr::filter(A==0 & delta1==1)
  new.data.A1T12 <- data.predictA1T12 <- data %>% dplyr::filter(A==1 & delta1==1)
  XmatA0 <- as.matrix(select(data.predictA0, Xnames))
  XmatA1 <- as.matrix(select(data.predictA1, Xnames))
  XmatA0T12 <- as.matrix(select(data.predictA0T12, Xnames))
  XmatA1T12 <- as.matrix(select(data.predictA1T12, Xnames))
  m.X.A0 <- apply(XmatA0, 2, mean)
  m.X.A1 <- apply(XmatA1, 2, mean)
  m.X.A0T12 <- apply(XmatA0T12, 2, mean)
  m.X.A1T12 <- apply(XmatA1T12, 2, mean)
  #### Formulas for the analysis
  if(Ltime==F)
  {
  formula01 <- as.formula(paste0("Surv(T1, delta1) ~", paste0(Xnames, collapse = "+"),
                                 " + offset(log.gamma)"))
  formula02 <- as.formula(paste0("Surv(T1, delta2only) ~", paste0(Xnames, collapse = "+"),
                                 " + offset(log.gamma)"))
  formula12 <- as.formula(paste0("Surv(T1, T2, delta2) ~",
                                 paste0(Xnames, collapse= "+") ,
                                 " + offset(log.gamma)"))
  }
  ########  Initial parameter values ########

  fit.a0.01 <- coxph(formula01, data = new.data.A0)
  fit.a0.02 <- coxph(formula02,  data = new.data.A0)
  fit.a0.12 <- coxph(formula12, data = new.data.A0T12)
  fit.a1.01 <- coxph(formula01, data = new.data.A1)
  fit.a1.02 <- coxph(formula02, data = new.data.A1)
  fit.a1.12 <- coxph(formula12, data = new.data.A1T12)
  est.beta.a0.01 <- coef(fit.a0.01)
  est.beta.a1.01 <- coef(fit.a1.01)
  est.beta.a0.02 <- coef(fit.a0.02)
  est.beta.a1.02 <- coef(fit.a1.02)
  est.beta.a0.12  <- coef(fit.a0.12)
  est.beta.a1.12  <- coef(fit.a1.12)
  old.betas <- new.betas <- c(est.beta.a0.01, est.beta.a1.01, est.beta.a0.02,
                              est.beta.a1.02, est.beta.a0.12, est.beta.a1.12)
  old.thetas <- new.thetas <- c(1,1)
  s.fit.a0.1 <- survfit(fit.a0.01, censor = FALSE)
  s.fit.a0.2 <- survfit(fit.a0.02, censor = FALSE)
  s.fit.a0.12 <- survfit(fit.a0.12, censor = FALSE)
  s.fit.a1.1 <- survfit(fit.a1.01, censor = FALSE)
  s.fit.a1.2 <- survfit(fit.a1.02,  censor = FALSE)
  s.fit.a1.12 <- survfit(fit.a1.12, censor = FALSE)
  step.A0T1 <- stepfun(x = s.fit.a0.1$time,
                       y = c(0, s.fit.a0.1$cumhaz*exp(-sum(est.beta.a0.01*m.X.A0)
                                                      -mean(new.data.A0$log.gamma))))
  step.A0T2 <- stepfun(x = s.fit.a0.2$time,
                       y = c(0, s.fit.a0.2$cumhaz*exp(- sum(est.beta.a0.02*m.X.A0)
                                                      - mean(new.data.A0$log.gamma))))
  step.A0T12 <- stepfun(x = s.fit.a0.12$time,
                        y = c(0, s.fit.a0.12$cumhaz*exp(- sum(est.beta.a0.12*m.X.A0T12)
                                                        - mean(new.data.A0T12$log.gamma))))
  step.A1T1 <- stepfun(x = s.fit.a1.1$time,
                       y = c(0, s.fit.a1.1$cumhaz*exp(-sum(est.beta.a1.01*m.X.A1)
                                                      -mean(new.data.A1$log.gamma))))
  step.A1T2 <- stepfun(x = s.fit.a1.2$time,
                       y = c(0, s.fit.a1.2$cumhaz*exp(- sum(est.beta.a1.02*m.X.A1)
                                                      - mean(new.data.A1$log.gamma))))
  step.A1T12 <- stepfun(x = s.fit.a1.12$time,
                        y = c(0, s.fit.a1.12$cumhaz*exp(-sum(est.beta.a1.12*m.X.A1T12)
                                                        - mean(new.data.A1T12$log.gamma))))
  ################################################################################
  ################################################################################################
  # ##### Calculate per-person cumulative hazards (with gamma=0)
  # ################################################################################################

  H.a0.01 <- step.A0T1(data.predictA0$T1) * exp(XmatA0%*%est.beta.a0.01)
  H.a0.02 <- step.A0T2(data.predictA0$T1) * exp(XmatA0%*%est.beta.a0.02)
  H.a0.12 <- step.A0T12(data.predictA0T12$T2) * exp(XmatA0T12%*%est.beta.a0.12)
  H.a0.12.T1 <- step.A0T12(data.predictA0T12$T1) * exp(XmatA0T12%*%est.beta.a0.12)
  H.a1.01 <- step.A1T1(data.predictA1$T1) * exp(XmatA1%*%est.beta.a1.01)
  H.a1.02 <- step.A1T2(data.predictA1$T1) * exp(XmatA1%*%est.beta.a1.02)
  H.a1.12 <- step.A1T12(data.predictA1T12$T2) * exp(XmatA1T12%*%est.beta.a1.12)
  H.a1.12.T1 <- step.A1T12(data.predictA1T12$T1) * exp(XmatA1T12%*%est.beta.a1.12)

  iter <- cond <- 0
  E.gamma <- E.log.gamma <- s.i <- vector(length = n.sample)
  ### Finally, the EM loop
  while(cond==0 & iter < max.iter)
  {
    iter <- iter + 1
    ##### E-step
    #####  Per-person posterior distriubtion parametrs
    s.i[A==0] <- H.a0.01 + H.a0.02
    s.i[A==1] <- H.a1.01 + H.a1.02
    s.i[A==0 & delta1==1] <- s.i[A==0 & delta1==1] + H.a0.12 - H.a0.12.T1
    s.i[A==1 & delta1==1] <- s.i[A==1 & delta1==1] + H.a1.12 - H.a1.12.T1
    E.gamma[A==0] <- (1/new.thetas[1] + delta1A0 + delta2A0) /
      (1/new.thetas[1] + s.i[A==0])
    E.gamma[A==1] <- (1/new.thetas[2] + delta1A1 + delta2A1) /
      (1/new.thetas[2] + s.i[A==1])
    E.log.gamma[A==0] <- digamma(1/new.thetas[1] + delta1A0 + delta2A0) -
      log(1/new.thetas[1] + s.i[A==0])
    E.log.gamma[A==1] <- digamma(1/new.thetas[2] + delta1A1 + delta2A1) -
      log(1/new.thetas[2] + s.i[A==1])
    data$log.gamma <- log(E.gamma)
    new.data.A0 <- data %>% dplyr::filter(A==0)
    new.data.A1 <- data %>% dplyr::filter(A==1)
    new.data.A0T12 <- data %>% dplyr::filter(A==0 & delta1==1)
    new.data.A1T12 <- data %>% dplyr::filter(A==1 & delta1==1)
    ###############################################################################################
    ##### M-step
    #### Conditonially on gamma, fit illness-death PH models
    ###############################################################################################
    fit.a0.01 <- coxph(formula01, data = new.data.A0)
    fit.a0.02 <- coxph(formula02,  data = new.data.A0)
    fit.a0.12 <- coxph(formula12, data = new.data.A0T12)
    fit.a1.01 <- coxph(formula01, data = new.data.A1)
    fit.a1.02 <- coxph(formula02, data = new.data.A1)
    fit.a1.12 <- coxph(formula12, data = new.data.A1T12)
    est.beta.a0.01 <- coef(fit.a0.01)
    est.beta.a1.01 <- coef(fit.a1.01)
    est.beta.a0.02 <- coef(fit.a0.02)
    est.beta.a1.02 <- coef(fit.a1.02)
    est.beta.a0.12  <- coef(fit.a0.12)
    est.beta.a1.12  <- coef(fit.a1.12)
    new.betas <- c(est.beta.a0.01, est.beta.a1.01, est.beta.a0.02,
                   est.beta.a1.02, est.beta.a0.12, est.beta.a1.12)
    ################################################################################################
    ##### Create step functions from all baseline hazard estimators ########
    ################################################################################################
    s.fit.a0.1 <- survfit(fit.a0.01, newdata = data.calc.H, censor = FALSE)
    s.fit.a0.2 <- survfit(fit.a0.02, newdata = data.calc.H, censor = FALSE)
    s.fit.a0.12 <- survfit(fit.a0.12, newdata = data.calc.H, censor = FALSE)
    s.fit.a1.1 <- survfit(fit.a1.01, newdata = data.calc.H, censor = FALSE)
    s.fit.a1.2 <- survfit(fit.a1.02, newdata = data.calc.H, censor = FALSE)
    s.fit.a1.12 <- survfit(fit.a1.12, newdata = data.calc.H, censor = FALSE)
    step.A0T1 <- stepfun(x = s.fit.a0.1$time, y = c(0, -log(s.fit.a0.1$surv)))
    step.A0T2 <- stepfun(x = s.fit.a0.2$time, y = c(0, -log(s.fit.a0.2$surv)))
    step.A0T12 <- stepfun(x = s.fit.a0.12$time, y = c(0, -log(s.fit.a0.12$surv)))
    step.A1T1 <- stepfun(x = s.fit.a1.1$time, y = c(0, -log(s.fit.a1.1$surv)))
    step.A1T2 <- stepfun(x = s.fit.a1.2$time, y = c(0, -log(s.fit.a1.2$surv)))
    step.A1T12 <- stepfun(x = s.fit.a1.12$time, y = c(0, -log(s.fit.a1.12$surv)))
    ################################################################################
    s.fit.a0.1 <- survfit(fit.a0.01, censor = FALSE)
    s.fit.a0.2 <- survfit(fit.a0.02, censor = FALSE)
    s.fit.a0.12 <- survfit(fit.a0.12, censor = FALSE)
    s.fit.a1.1 <- survfit(fit.a1.01, censor = FALSE)
    s.fit.a1.2 <- survfit(fit.a1.02,  censor = FALSE)
    s.fit.a1.12 <- survfit(fit.a1.12, censor = FALSE)
    step.A0T1 <- stepfun(x = s.fit.a0.1$time,
                         y = c(0, s.fit.a0.1$cumhaz*exp(-sum(est.beta.a0.01*m.X.A0)
                                                        -mean(new.data.A0$log.gamma))))
    step.A0T2 <- stepfun(x = s.fit.a0.2$time,
                         y = c(0, s.fit.a0.2$cumhaz*exp(- sum(est.beta.a0.02*m.X.A0)
                                                        - mean(new.data.A0$log.gamma))))

    step.A0T12 <- stepfun(x = s.fit.a0.12$time,
                          y = c(0, s.fit.a0.12$cumhaz*exp(- sum(est.beta.a0.12*m.X.A0T12)
                                                          - mean(new.data.A0T12$log.gamma))))
    step.A1T1 <- stepfun(x = s.fit.a1.1$time,
                         y = c(0, s.fit.a1.1$cumhaz*exp(-sum(est.beta.a1.01*m.X.A1)
                                                        -mean(new.data.A1$log.gamma))))
    step.A1T2 <- stepfun(x = s.fit.a1.2$time,
                         y = c(0, s.fit.a1.2$cumhaz*exp(- sum(est.beta.a1.02*m.X.A1)
                                                        - mean(new.data.A1$log.gamma))))
    step.A1T12 <- stepfun(x = s.fit.a1.12$time,
                          y = c(0, s.fit.a1.12$cumhaz*exp(-sum(est.beta.a1.12*m.X.A1T12)
                                                          - mean(new.data.A1T12$log.gamma))))
    # Per-person cumulative hazards (with gamma=0) for "posterior distriubtion"
    ##################################################################################################
    H.a0.01 <- step.A0T1(data.predictA0$T1) * exp(XmatA0%*%est.beta.a0.01)
    H.a0.02 <- step.A0T2(data.predictA0$T1) * exp(XmatA0%*%est.beta.a0.02)
    H.a0.12 <- step.A0T12(data.predictA0T12$T2) * exp(XmatA0T12%*%est.beta.a0.12)
    H.a0.12.T1 <- step.A0T12(data.predictA0T12$T1) * exp(XmatA0T12%*%est.beta.a0.12)
    H.a1.01 <- step.A1T1(data.predictA1$T1) * exp(XmatA1%*%est.beta.a1.01)
    H.a1.02 <- step.A1T2(data.predictA1$T1) * exp(XmatA1%*%est.beta.a1.02)
    H.a1.12 <- step.A1T12(data.predictA1T12$T2) * exp(XmatA1T12%*%est.beta.a1.12)
    H.a1.12.T1 <- step.A1T12(data.predictA1T12$T1) * exp(XmatA1T12%*%est.beta.a1.12)
    ##### New theta values #####
    new.thetas[1] <- optimize(f = Eloglik, interval = c(0.01, 30),
                              delta1 = delta1A0, delta2 = delta2A0,
                              E.gamma = E.gamma[A==0],
                              E.log.gamma = E.log.gamma[A==0],
                              maximum = T)$maximum
    new.thetas[2] <- optimize(f = Eloglik, interval = c(0.01, 30),
                              delta1 = delta1A1, delta2 = delta2A1,
                              E.gamma = E.gamma[A==1],
                              E.log.gamma = E.log.gamma[A==1],
                              maximum = T)$maximum
    if (max(abs(c(new.thetas - old.thetas, new.betas - old.betas))) < 0.0001)
    {
      cond <- 1
    }
    old.betas <- new.betas
    old.thetas <- new.thetas
  }
  fit.list <- list(fit.a0.01 = fit.a0.01, fit.a0.02 = fit.a0.02, fit.a0.12 = fit.a0.12,
                   fit.a1.01 = fit.a1.01, fit.a1.02 = fit.a1.02, fit.a1.12 = fit.a1.12)
  list.out <- list(betas = new.betas, thetas = new.thetas,
                   fit.list = fit.list, iter = iter)
  return(list.out)
}
