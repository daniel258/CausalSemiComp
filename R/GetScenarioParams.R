# A function to pick scenarios for simulations
# Scenario kxxxy means we use "scenario k" with correlation xxx/100 and Kendall's tau 0.25 (y=1), 0.5 (y=2) or 0.75 (y=3)
# Scenarios 1xxxy-3xxxy for the nonparametric bounds.
# Scenarios 1xxxy for adjusted bounds and beta values
# Scenario 2xxxy and 3xxxy for null effect, but different strata proportions

GetScenarioParams <- function(scenario.num)
{
  if (floor(scenario.num/10000)==1)
  {
    base.weib.scale.a0.01 <- 12.5
    base.weib.scale.a1.01 <- 10
    base.weib.scale.a0.02 <- 25
    base.weib.scale.a1.02 <- 17.5
    base.weib.scale.a0.12 <- 25
    base.weib.scale.a1.12 <- 20
    base.weib.shape.a0.01 <- 2
    base.weib.shape.a1.01 <- 3
    base.weib.shape.a0.02 <- 2.25
    base.weib.shape.a1.02 <- 1.5
    base.weib.shape.a0.12 <- 1.5
    base.weib.shape.a1.12 <- 2.5
    beta.a0.01 <- log(c(0.25, 1))
    beta.a0.02 <- log(c(1, 1))
    beta.a0.12 <- log(c(1, 1))
    beta.a1.01 <- log(c(1, 1))
    beta.a1.02 <- log(c(1, 1))
    beta.a1.12 <- log(c(1, 1))
    tau <- (scenario.num %% 10) * 0.25
    theta <- 2*tau/(1 - tau)
    rho <- round(scenario.num/10)/100 - 10
  }
  if (floor(scenario.num/10000)==2)
  {
    base.weib.scale.a0.01 <- 20
    base.weib.scale.a1.01 <- 20
    base.weib.scale.a0.02 <- 15
    base.weib.scale.a1.02 <- 15
    base.weib.scale.a0.12 <- 25
    base.weib.scale.a1.12 <- 25
    base.weib.shape.a0.01 <- 3
    base.weib.shape.a1.01 <- 3
    base.weib.shape.a0.02 <- 1.5
    base.weib.shape.a1.02 <- 1.5
    base.weib.shape.a0.12 <- 2.75
    base.weib.shape.a1.12 <- 2.75
    beta.a0.01 <- log(c(1, 1))
    beta.a0.02 <- log(c(1, 1))
    beta.a0.12 <- log(c(1, 1))
    beta.a1.01 <- log(c(1, 1))
    beta.a1.02 <- log(c(1, 1))
    beta.a1.12 <- log(c(1, 1))
    beta.a0.01 <- log(c(1, 1))
    beta.a0.02 <- log(c(1, 1))
    beta.a0.12 <- log(c(1, 1))
    beta.a1.01 <- log(c(1, 1))
    beta.a1.02 <- log(c(1, 1))
    beta.a1.12 <- log(c(1, 1))
    tau <- (scenario.num %% 10) * 0.25
    theta <- 2*tau/(1 - tau)
    rho <- round(scenario.num/10)/100 - 20
  }
  if (floor(scenario.num/10000)==3)
  {
    base.weib.scale.a0.01 <- 7.5
    base.weib.scale.a1.01 <- 7.5
    base.weib.scale.a0.02 <- 15
    base.weib.scale.a1.02 <- 15
    base.weib.scale.a0.12 <- 20
    base.weib.scale.a1.12 <- 20
    base.weib.shape.a0.01 <- 2
    base.weib.shape.a1.01 <- 2
    base.weib.shape.a0.02 <- 1.75
    base.weib.shape.a1.02 <- 1.75
    base.weib.shape.a0.12 <- 2.5
    base.weib.shape.a1.12 <- 2.5
    beta.a0.01 <- log(c(1, 1))
    beta.a0.02 <- log(c(1, 1))
    beta.a0.12 <- log(c(1, 1))
    beta.a1.01 <- log(c(1, 1))
    beta.a1.02 <- log(c(1, 1))
    beta.a1.12 <- log(c(1, 1))
    tau <- (scenario.num %% 10) * 0.25
    theta <- 2*tau/(1 - tau)
    rho <- round(scenario.num/10)/100 - 30
  }

  if (floor(scenario.num/10000)==4)
  {
    base.weib.scale.a0.01 <- 5
    base.weib.scale.a1.01 <- 2.5
    base.weib.scale.a0.02 <- 7.5
    base.weib.scale.a1.02 <- 5
    base.weib.scale.a0.12 <- 15
    base.weib.scale.a1.12 <- 10
    base.weib.shape.a0.01 <- 2.5
    base.weib.shape.a1.01 <- 2.5
    base.weib.shape.a0.02 <- 1.5
    base.weib.shape.a1.02 <- 1.5
    base.weib.shape.a0.12 <- 2.5
    base.weib.shape.a1.12 <- 2.5
    beta.a0.01 <- log(c(0.25, 1.5))
    beta.a0.02 <- log(c(0.75, 1.5))
    beta.a0.12 <- log(c(1, 1))
    beta.a1.01 <- log(c(1, 1.5))
    beta.a1.02 <- log(c(1, 2.5))
    beta.a1.12 <- log(c(0.5, 2))
    tau <- (scenario.num %% 10) * 0.25
    theta <- 2*tau/(1 - tau)
    rho <- round(scenario.num/10)/100 - 40
  }
  if (floor(scenario.num/10000)==5)
  {
    base.weib.scale.a0.01 <- 10
    base.weib.scale.a1.01 <- 10
    base.weib.scale.a0.02 <- 15
    base.weib.scale.a1.02 <- 15
    base.weib.scale.a0.12 <- 25
    base.weib.scale.a1.12 <- 25
    base.weib.shape.a0.01 <- 3
    base.weib.shape.a1.01 <- 3
    base.weib.shape.a0.02 <- 1.5
    base.weib.shape.a1.02 <- 1.5
    base.weib.shape.a0.12 <- 2.75
    base.weib.shape.a1.12 <- 2.75
    beta.a0.01 <- log(c(1, 1))
    beta.a0.02 <- log(c(1, 1))
    beta.a0.12 <- log(c(1, 1))
    beta.a1.01 <- log(c(1, 1))
    beta.a1.02 <- log(c(1, 1))
    beta.a1.12 <- log(c(1, 1))
    tau <- (scenario.num %% 10) * 0.25
    theta <- 2*tau/(1 - tau)
    rho <- round(scenario.num/10)/100 - 50
  }
  params <- mget(ls(environment()))
  return(params)
}
