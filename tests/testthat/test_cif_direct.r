library(survival)
library(riskRegression)
library(prodlim)

set.seed(42)

sim_dat <- adjustedCurves::sim_confounded_surv(n=50)
sim_dat$event[sim_dat$event==1] <- sample(c(1, 2), size=sum(sim_dat$event),
                                          replace=T)
sim_dat$group <- as.factor(sim_dat$group)

# outcome model
mod <- riskRegression::CSC(Hist(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6,
                           data=sim_dat)
mod2 <- riskRegression::FGR(Hist(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6,
                            data=sim_dat, cause=1)

## Just check if function throws any errors
test_that("CSC, 2 treatments, no conf_int, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=FALSE,
                                           outcome_model=mod,
                                           cause=1), NA)
})

test_that("CSC, 2 treatments, with conf_int, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=TRUE,
                                           outcome_model=mod,
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
                                           outcome_model=mod,
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
                                           outcome_model=mod,
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
                                           outcome_model=mod,
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
                                           outcome_model=mod2,
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
                                           outcome_model=mod2,
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
                                           outcome_model=mod2,
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
                                           outcome_model=mod2,
                                           times=c(0.5, 1, 1.2),
                                           cause=1), NA)
})

set.seed(34)
sim_dat <- adjustedCurves::sim_confounded_surv(n=120)
sim_dat$group[sim_dat$group==1] <- sample(c(1, 2),
                                          size=nrow(sim_dat[sim_dat$group==1,]),
                                          replace=T)
sim_dat$group <- as.factor(sim_dat$group)
sim_dat$event[sim_dat$event==1] <- sample(c(1, 2), size=sum(sim_dat$event),
                                          replace=T)

# outcome models
mod <- riskRegression::CSC(Hist(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6,
                           data=sim_dat)
mod2 <- riskRegression::FGR(Hist(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6,
                            data=sim_dat, cause=1)


test_that("CSC, > 2 treatments, no conf_int, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=FALSE,
                                           outcome_model=mod,
                                           cause=1), NA)
})

test_that("CSC, > 2 treatments, with conf_int, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=TRUE,
                                           outcome_model=mod,
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
                                           outcome_model=mod,
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
                                           outcome_model=mod,
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
                                           outcome_model=mod,
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
                                           outcome_model=mod,
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
                                           outcome_model=mod2,
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
                                           outcome_model=mod2,
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
                                           outcome_model=mod2,
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
                                           outcome_model=mod2,
                                           times=c(0.5, 1, 1.2),
                                           cause=1), NA)
})
