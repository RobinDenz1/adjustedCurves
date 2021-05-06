library(adjustedCurves)
library(survival)

set.seed(42)

sim_dat <- adjustedCurves::sim_confounded_surv(n=100)
sim_dat$group <- as.factor(sim_dat$group)

outcome_vars <- c("x1", "x2", "x3", "x4")

## Geese
test_that("2 treatments, no conf_int, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            model_type="geese",
                                            method="direct_pseudo",
                                            conf_int=F,
                                            outcome_vars=outcome_vars,
                                            type_time="bs"), NA)
})

test_that("2 treatments, no conf_int, with boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            model_type="geese",
                                            method="direct_pseudo",
                                            conf_int=F,
                                            bootstrap=T,
                                            n_boot=2,
                                            outcome_vars=outcome_vars,
                                            type_time="bs"), NA)
})

test_that("2 treatments, no conf_int, no boot, with times", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            model_type="geese",
                                            method="direct_pseudo",
                                            conf_int=F,
                                            bootstrap=T,
                                            n_boot=2,
                                            outcome_vars=outcome_vars,
                                            times=c(0.3, 0.8),
                                            type_time="factor"), NA)
})

test_that("2 treatments, no conf_int, no boot, with times, type_time='factor'", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            model_type="geese",
                                            method="direct_pseudo",
                                            conf_int=F,
                                            bootstrap=F,
                                            n_boot=5,
                                            outcome_vars=outcome_vars,
                                            times=c(0.3, 0.8),
                                            type_time="factor"), NA)
})

## series of linear models
test_that("2 treatments, no conf_int, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            model_type="lm",
                                            method="direct_pseudo",
                                            conf_int=F,
                                            outcome_vars=outcome_vars), NA)
})

test_that("2 treatments, no conf_int, with boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            model_type="lm",
                                            method="direct_pseudo",
                                            conf_int=F,
                                            bootstrap=T,
                                            n_boot=2,
                                            outcome_vars=outcome_vars), NA)
})

test_that("2 treatments, no conf_int, no boot, with times", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            model_type="lm",
                                            method="direct_pseudo",
                                            conf_int=F,
                                            bootstrap=T,
                                            n_boot=2,
                                            outcome_vars=outcome_vars,
                                            times=c(0.3, 0.8)), NA)
})

test_that("2 treatments, conf_int, no boot, with times", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            model_type="lm",
                                            method="direct_pseudo",
                                            conf_int=T,
                                            bootstrap=F,
                                            n_boot=5,
                                            outcome_vars=outcome_vars,
                                            times=c(0.3, 0.8),
                                            type_time="factor"), NA)
})

sim_dat <- adjustedCurves::sim_confounded_surv(n=100)
sim_dat$group[sim_dat$group==1] <- sample(c(1, 2),
                                          size=nrow(sim_dat[sim_dat$group==1,]),
                                          replace=T)
sim_dat$group <- as.factor(sim_dat$group)
sim_dat$time <- round(sim_dat$time, 1)

## geese
test_that("> 2 treatments, no conf_int, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            model_type="geese",
                                            method="direct_pseudo",
                                            conf_int=F,
                                            outcome_vars=outcome_vars,
                                            type_time="bs"), NA)
})

test_that("> 2 treatments, no conf_int, with boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            model_type="geese",
                                            method="direct_pseudo",
                                            conf_int=F,
                                            bootstrap=T,
                                            n_boot=2,
                                            outcome_vars=outcome_vars,
                                            type_time="bs"), NA)
})

test_that("> 2 treatments, no conf_int, no boot, with times", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            model_type="geese",
                                            method="direct_pseudo",
                                            conf_int=F,
                                            bootstrap=F,
                                            n_boot=2,
                                            outcome_vars=outcome_vars,
                                            times=c(0.3, 0.8, 1, 1.5),
                                            type_time="factor"), NA)
})

test_that("> 2 treatments, no conf_int, no boot, with times, type_time='factor'", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            model_type="geese",
                                            method="direct_pseudo",
                                            conf_int=F,
                                            bootstrap=F,
                                            n_boot=5,
                                            outcome_vars=outcome_vars,
                                            times=c(0.3, 0.8, 1, 1.5),
                                            type_time="factor"), NA)
})

## series of linear models
test_that("2 treatments, no conf_int, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            model_type="lm",
                                            method="direct_pseudo",
                                            conf_int=F,
                                            outcome_vars=outcome_vars), NA)
})

test_that("2 treatments, no conf_int, with boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            model_type="lm",
                                            method="direct_pseudo",
                                            conf_int=F,
                                            bootstrap=T,
                                            n_boot=2,
                                            outcome_vars=outcome_vars), NA)
})

test_that("2 treatments, no conf_int, no boot, with times", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            model_type="lm",
                                            method="direct_pseudo",
                                            conf_int=F,
                                            bootstrap=T,
                                            n_boot=2,
                                            outcome_vars=outcome_vars,
                                            times=c(0.3, 0.8)), NA)
})

test_that("2 treatments, conf_int, no boot, with times", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            model_type="lm",
                                            method="direct_pseudo",
                                            conf_int=T,
                                            bootstrap=F,
                                            n_boot=5,
                                            outcome_vars=outcome_vars,
                                            times=c(0.3, 0.8),
                                            type_time="factor"), NA)
})
