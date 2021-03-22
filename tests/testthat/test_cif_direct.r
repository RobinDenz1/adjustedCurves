library(survival)
library(riskRegression)

set.seed(42)

sim_dat <- adjustedCurves::sim_confounded_surv(n=100)
sim_dat$event[sim_dat$event==1] <- sample(c(1, 2), size=sum(sim_dat$event),
                                          replace=T)
sim_dat$group <- as.factor(sim_dat$group)

# outcome model
mod <- riskRegression::CSC(Hist(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6,
                           data=sim_dat)

## Just check if function throws any errors
test_that("2 treatments, no conf_int, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=F,
                                           outcome_model=mod,
                                           cause=1), NA)
})

test_that("2 treatments, with conf_int, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=T,
                                           outcome_model=mod,
                                           cause=1), NA)
})

test_that("2 treatments, no conf_int, with boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=F,
                                           bootstrap=T,
                                           n_boot=2,
                                           outcome_model=mod,
                                           cause=1), NA)
})

test_that("2 treatments, with conf_int, with boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=T,
                                           bootstrap=T,
                                           n_boot=2,
                                           outcome_model=mod,
                                           cause=1), NA)
})

test_that("2 treatments, no conf_int, no boot, with times", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=T,
                                           bootstrap=T,
                                           n_boot=10,
                                           outcome_model=mod,
                                           times=c(0.5, 1, 1.2),
                                           cause=1), NA)
})

sim_dat <- adjustedCurves::sim_confounded_surv(n=300)
sim_dat$group[sim_dat$group==1] <- sample(c(1, 2),
                                          size=nrow(sim_dat[sim_dat$group==1,]),
                                          replace=T)
sim_dat$group <- as.factor(sim_dat$group)
sim_dat$event[sim_dat$event==1] <- sample(c(1, 2), size=sum(sim_dat$event),
                                          replace=T)

# outcome model
mod <- riskRegression::CSC(Hist(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6,
                           data=sim_dat)


test_that("> 2 treatments, no conf_int, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=F,
                                           outcome_model=mod,
                                           cause=1), NA)
})

test_that("> 2 treatments, with conf_int, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=T,
                                           outcome_model=mod,
                                           cause=1), NA)
})

test_that("> 2 treatments, no conf_int, with boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=F,
                                           bootstrap=T,
                                           n_boot=2,
                                           outcome_model=mod,
                                           cause=1), NA)
})

test_that("> 2 treatments, with conf_int, with boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=T,
                                           bootstrap=T,
                                           n_boot=2,
                                           outcome_model=mod,
                                           cause=1), NA)
})

test_that("> 2 treatments, no conf_int, no boot, with times", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=F,
                                           bootstrap=F,
                                           n_boot=2,
                                           outcome_model=mod,
                                           times=c(0.5, 1, 1.2),
                                           cause=1), NA)
})



