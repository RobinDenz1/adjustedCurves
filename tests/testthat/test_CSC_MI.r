library(survival)
library(riskRegression)
library(mice)

set.seed(42)

sim_dat <- adjustedCurves::sim_confounded_crisk(n=100)
sim_dat$group <- as.factor(sim_dat$group)

sim_dat$x1 <- ifelse(stats::runif(n=100)<0.5, NA, sim_dat$x1)

mids <- mice::mice(data=sim_dat, m=3, method="pmm", printFlag=F)

## Just check if function throws any errors
test_that("no additional args", {
  expect_error(adjustedCurves::CSC_MI(mids=mids,
                                      formula=Hist(time, event) ~ x1 + x2), NA)
})


