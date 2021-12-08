library(survival)
library(riskRegression)
library(prodlim)
library(mice)

set.seed(42)

sim_dat <- adjustedCurves::sim_confounded_crisk(n=100)
sim_dat$group <- as.factor(sim_dat$group)

sim_dat$x1 <- ifelse(stats::runif(n=100)<0.5, NA, sim_dat$x1)

mids <- mice::mice(data=sim_dat, m=3, method="pmm", printFlag=FALSE)

test_that("no additional args", {
  expect_error(adjustedCurves::FGR_MI(mids=mids,
                                      formula=Hist(time, event) ~ x1 + x2), NA)
})

test_that("with additional args", {
  expect_error(adjustedCurves::FGR_MI(mids=mids,
                                      formula=Hist(time, event) ~ x1 + x2,
                                      cause=1), NA)
})

test_that("trying to use data", {
  expect_error(adjustedCurves::FGR_MI(mids=mids,
                                      formula=Hist(time, event) ~ x1 + x2,
                                      cause=1,
                                      data=sim_dat),
               NULL)
})
