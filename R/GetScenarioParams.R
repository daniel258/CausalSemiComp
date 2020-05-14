# A function to pick scenarios for simulations
# Scenarios 1-3 for the nonparametric bounds.
# Scenario 1 for adjusted bounds and beta values
# Senarios 2-3 for
GetScenarioParams <- function(scenario.num)
{
  if (scenario.num==1)
  {
    base.weib.scale.a0.01 <- 12.5
    base.weib.scale.a1.01 <- 7.5
    base.weib.scale.a0.02 <- 20
    base.weib.scale.a1.02 <- 15
    base.weib.scale.a0.12 <- 25
    base.weib.scale.a1.12 <- 15
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
    theta <- 6
    rho <- 0.7
  }
  if (scenario.num==2)
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
    theta <- 2
  }
  if (scenario.num==3)
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
    theta <- 2
  }
  if (scenario.num==4)
  {
    base.weib.scale.a0.01 <- 10
    base.weib.scale.a1.01 <- 5
    base.weib.scale.a0.02 <- 20
    base.weib.scale.a1.02 <- 15
    base.weib.scale.a0.12 <- 25
    base.weib.scale.a1.12 <- 20
    base.weib.shape.a0.01 <- 2.75
    base.weib.shape.a1.01 <- 2.75
    base.weib.shape.a0.02 <- 1.5
    base.weib.shape.a1.02 <- 1.5
    base.weib.shape.a0.12 <- 2.5
    base.weib.shape.a1.12 <- 2.5
    beta.a0.01 <- log(c(0.3, 1))
    beta.a0.02 <- log(c(1, 1))
    beta.a0.12 <- log(c(1, 1))
    beta.a1.01 <- log(c(1, 1))
    beta.a1.02 <- log(c(0.5, 1))
    beta.a1.12 <- log(c(0.5, 1))
    theta <- 4
  }
  if (scenario.num==5)
  {
    base.weib.scale.a0.01 <- 12.5
    base.weib.scale.a1.01 <- 7.5
    base.weib.scale.a0.02 <- 15
    base.weib.scale.a1.02 <- 15
    base.weib.scale.a0.12 <- 20
    base.weib.scale.a1.12 <- 10
    base.weib.shape.a0.01 <- 2.5
    base.weib.shape.a1.01 <- 2.5
    base.weib.shape.a0.02 <- 2
    base.weib.shape.a1.02 <- 1.5
    base.weib.shape.a0.12 <- 1.5
    base.weib.shape.a1.12 <- 2
    beta.a0.01 <- log(c(1, 1))
    beta.a0.02 <- log(c(1, 1))
    beta.a0.12 <- log(c(1, 1))
    beta.a1.01 <- log(c(1, 1))
    beta.a1.02 <- log(c(1, 1))
    beta.a1.12 <- log(c(1, 1))
    theta <- 2
  }
  if (scenario.num==6)
  {
    base.weib.scale.a0.01 <- 10
    base.weib.scale.a1.01 <- 7.5
    base.weib.scale.a0.02 <- 20
    base.weib.scale.a1.02 <- 17.5
    base.weib.scale.a0.12 <- 30
    base.weib.scale.a1.12 <- 17.5
    base.weib.shape.a0.01 <- 2.5
    base.weib.shape.a1.01 <- 3
    base.weib.shape.a0.02 <- 1.5
    base.weib.shape.a1.02 <- 1.5
    base.weib.shape.a0.12 <- 2.5
    base.weib.shape.a1.12 <- 2.5
    beta.a0.01 <- log(c(0.3, 1))
    beta.a0.02 <- log(c(1, 1))
    beta.a0.12 <- log(c(1, 1))
    beta.a1.01 <- log(c(0.8, 1))
    beta.a1.02 <- log(c(0.8, 1))
    beta.a1.12 <- log(c(0.5, 1))
    theta <- 6
  }
  if (scenario.num==7)
  {
    base.weib.scale.a0.01 <- 10
    base.weib.scale.a1.01 <- 7.5
    base.weib.scale.a0.02 <- 20
    base.weib.scale.a1.02 <- 17.5
    base.weib.scale.a0.12 <- 30
    base.weib.scale.a1.12 <- 17.5
    base.weib.shape.a0.01 <- 1
    base.weib.shape.a1.01 <- 1
    base.weib.shape.a0.02 <- 1
    base.weib.shape.a1.02 <- 1
    base.weib.shape.a0.12 <- 1
    base.weib.shape.a1.12 <- 1
    beta.a0.01 <- log(c(0.3, 1))
    beta.a0.02 <- log(c(1, 1))
    beta.a0.12 <- log(c(1, 1))
    beta.a1.01 <- log(c(0.8, 1))
    beta.a1.02 <- log(c(0.8, 1))
    beta.a1.12 <- log(c(0.5, 1))
    theta <- 2
  }
  params <- mget(ls(environment()))
  return(params)
}
