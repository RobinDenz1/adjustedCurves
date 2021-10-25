library(survival)
library(riskRegression)
library(prodlim)

set.seed(42)

sim_dat <- adjustedCurves::sim_confounded_surv(n=65)
sim_dat$event[sim_dat$event==1] <- sample(c(1, 2), size=sum(sim_dat$event),
                                          replace=TRUE)
sim_dat$group <- as.factor(sim_dat$group)

# outcome model
mod_CSC <- riskRegression::CSC(Hist(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6,
                           data=sim_dat)
mod_FGR <- riskRegression::FGR(Hist(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6,
                            data=sim_dat, cause=1)
# fit models
mod_riskRegression <- riskRegression::riskRegression(Hist(time, event) ~ x1 +
                                                       x2 + x3 + group,
                                                     data=sim_dat, cause=1)
mod_ARR <- riskRegression::ARR(Hist(time, event) ~ x1 + x2 + x3 + group,
                               data=sim_dat, cause=1)
mod_prodlim <- prodlim::prodlim(prodlim::Hist(time, event) ~ x1 + x2 + x3 +
                                  group, data=sim_dat)

## Just check if function throws any errors
test_that("CSC, 2 treatments, no conf_int, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=FALSE,
                                           outcome_model=mod_CSC,
                                           cause=1), NA)
})

test_that("CSC, 2 treatments, with conf_int, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=TRUE,
                                           outcome_model=mod_CSC,
                                           cause=1), NA)
})

test_that("CSC, 2 treatments, no conf_int, with boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=FALSE,
                                           bootstrap=TRUE,
                                           n_boot=2,
                                           outcome_model=mod_CSC,
                                           cause=1), NA)
})

test_that("CSC, 2 treatments, with conf_int, with boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=TRUE,
                                           bootstrap=TRUE,
                                           n_boot=2,
                                           outcome_model=mod_CSC,
                                           cause=1), NA)
})

test_that("CSC, 2 treatments, no conf_int, no boot, with times", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=TRUE,
                                           bootstrap=TRUE,
                                           n_boot=2,
                                           outcome_model=mod_CSC,
                                           times=c(0.5, 1, 1.2),
                                           cause=1), NA)
})

test_that("FGR, 2 treatments, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=FALSE,
                                           outcome_model=mod_FGR,
                                           cause=1), NA)
})

test_that("FGR, 2 treatments, with boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=FALSE,
                                           bootstrap=TRUE,
                                           n_boot=2,
                                           outcome_model=mod_FGR,
                                           cause=1), NA)
})

test_that("FGR, 2 treatments, no boot, with times", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=FALSE,
                                           bootstrap=FALSE,
                                           n_boot=2,
                                           outcome_model=mod_FGR,
                                           times=c(0.5, 1, 1.2),
                                           cause=1), NA)
})

test_that("FGR, 2 treatments, with boot, with times", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=FALSE,
                                           bootstrap=TRUE,
                                           n_boot=2,
                                           outcome_model=mod_FGR,
                                           times=c(0.5, 1, 1.2),
                                           cause=1), NA)
})

test_that("riskRegression, 2 treatments, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=FALSE,
                                           outcome_model=mod_riskRegression,
                                           cause=1), NA)
})

test_that("riskRegression, 2 treatments, with boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=TRUE,
                                           outcome_model=mod_riskRegression,
                                           cause=1), NA)
})

test_that("ARR, 2 treatments, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=FALSE,
                                           outcome_model=mod_ARR,
                                           cause=1), NA)
})

test_that("ARR, 2 treatments, with boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=TRUE,
                                           outcome_model=mod_ARR,
                                           cause=1), NA)
})

test_that("prodlim, 2 treatments, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=FALSE,
                                           outcome_model=mod_prodlim,
                                           cause=1), NA)
})

test_that("prodlim, 2 treatments, with boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=TRUE,
                                           outcome_model=mod_prodlim,
                                           cause=1), NA)
})

set.seed(34)
sim_dat <- adjustedCurves::sim_confounded_surv(n=160)
sim_dat$group[sim_dat$group==1] <- sample(c(1, 2),
                                          size=nrow(sim_dat[sim_dat$group==1,]),
                                          replace=T)
sim_dat$group <- as.factor(sim_dat$group)
sim_dat$event[sim_dat$event==1] <- sample(c(1, 2), size=sum(sim_dat$event),
                                          replace=TRUE)

# outcome models
mod_CSC <- riskRegression::CSC(Hist(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6,
                           data=sim_dat)
mod_FGR <- riskRegression::FGR(Hist(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6,
                            data=sim_dat, cause=1)
mod_riskRegression <- riskRegression::riskRegression(Hist(time, event) ~ x1 +
                                                       x2 + x3 + group,
                                                     data=sim_dat, cause=1)
mod_ARR <- riskRegression::ARR(Hist(time, event) ~ x1 + x2 + x3 + group,
                               data=sim_dat, cause=1)
mod_prodlim <- prodlim::prodlim(prodlim::Hist(time, event) ~ x1 + x2 + x3 +
                                  group, data=sim_dat)


test_that("CSC, > 2 treatments, no conf_int, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=FALSE,
                                           outcome_model=mod_CSC,
                                           cause=1), NA)
})

test_that("CSC, > 2 treatments, with conf_int, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=TRUE,
                                           outcome_model=mod_CSC,
                                           cause=1), NA)
})

test_that("CSC, > 2 treatments, no conf_int, with boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=FALSE,
                                           bootstrap=TRUE,
                                           n_boot=2,
                                           outcome_model=mod_CSC,
                                           cause=1), NA)
})

test_that("CSC, > 2 treatments, with conf_int, with boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=TRUE,
                                           bootstrap=TRUE,
                                           n_boot=2,
                                           outcome_model=mod_CSC,
                                           cause=1), NA)
})

test_that("CSC, > 2 treatments, no conf_int, no boot, with times", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=FALSE,
                                           bootstrap=FALSE,
                                           n_boot=2,
                                           outcome_model=mod_CSC,
                                           times=c(0.5, 1, 1.2),
                                           cause=1), NA)
})

test_that("CSC, > 2 treatments, no conf_int, with boot, with times", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=FALSE,
                                           bootstrap=TRUE,
                                           n_boot=2,
                                           outcome_model=mod_CSC,
                                           times=c(0.5, 1, 1.2),
                                           cause=1), NA)
})

test_that("FGR, > 2 treatments, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=FALSE,
                                           outcome_model=mod_FGR,
                                           cause=1), NA)
})

test_that("FGR, > 2 treatments, with boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=FALSE,
                                           bootstrap=TRUE,
                                           n_boot=2,
                                           outcome_model=mod_FGR,
                                           cause=1), NA)
})

test_that("FGR, > 2 treatments, no boot, with times", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=FALSE,
                                           bootstrap=FALSE,
                                           n_boot=2,
                                           outcome_model=mod_FGR,
                                           times=c(0.5, 1, 1.2),
                                           cause=1), NA)
})

test_that("FGR, > 2 treatments, with boot, with times", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=FALSE,
                                           bootstrap=TRUE,
                                           n_boot=2,
                                           outcome_model=mod_FGR,
                                           times=c(0.5, 1, 1.2),
                                           cause=1), NA)
})

test_that("riskRegression, > 2 treatments, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=FALSE,
                                           outcome_model=mod_riskRegression,
                                           cause=1), NA)
})

test_that("riskRegression, > 2 treatments, with boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=TRUE,
                                           outcome_model=mod_riskRegression,
                                           cause=1), NA)
})

test_that("ARR, > 2 treatments, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=FALSE,
                                           outcome_model=mod_ARR,
                                           cause=1), NA)
})

test_that("ARR, > 2 treatments, with boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=TRUE,
                                           outcome_model=mod_ARR,
                                           cause=1), NA)
})

test_that("prodlim, > 2 treatments, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=FALSE,
                                           outcome_model=mod_prodlim,
                                           cause=1), NA)
})

test_that("prodlim, > 2 treatments, with boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=TRUE,
                                           outcome_model=mod_prodlim,
                                           cause=1), NA)
})

# NOTE: These models are supported, but would require all kinds of package
#       dependencies if tested here. The tests are run in the development
#       process but are commented out here.

#library(randomForestSRC)
#mod_rfsrc <- randomForestSRC::rfsrc(Surv(time, event) ~ x1 + x2 + x3 + group,
#                                   data=sim_dat)
#library(casebase)
#mod_fitSmoothHazard <- casebase::fitSmoothHazard(event ~ time + x1 + group,
#                                                 sim_dat, ratio=10)
#
#
#test_that("rfsrc, 2 treatments, no boot", {
#  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
#                                           variable="group",
#                                           ev_time="time",
#                                           event="event",
#                                           method="direct",
#                                           conf_int=FALSE,
#                                           outcome_model=mod_rfsrc,
#                                           cause=1), NA)
#})
#
#test_that("rfsrc, 2 treatments, with boot", {
#  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
#                                           variable="group",
#                                           ev_time="time",
#                                           event="event",
#                                           method="direct",
#                                           conf_int=TRUE,
#                                           outcome_model=mod_rfsrc,
#                                           cause=1), NA)
#})
#
## NOTE: Currently doesn't work due to bugs in casebase
#test_that("fitSmoothHazard, 2 treatments, no boot", {
#  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
#                                           variable="group",
#                                           ev_time="time",
#                                           event="event",
#                                           method="direct",
#                                           conf_int=FALSE,
#                                           outcome_model=mod_fitSmoothHazard,
#                                           cause=1), NA)
#})
#
## NOTE Currently doesn't work due to bugs in casebase
#test_that("fitSmoothHazard, 2 treatments, with boot", {
#  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
#                                           variable="group",
#                                           ev_time="time",
#                                           event="event",
#                                           method="direct",
#                                           conf_int=TRUE,
#                                           outcome_model=mod_fitSmoothHazard,
#                                           cause=1), NA)
#})
