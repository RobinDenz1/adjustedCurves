#library(adjustedCurves)
library(survival)

set.seed(42)

sim_dat <- adjustedCurves::sim_confounded_surv(n=200)
sim_dat$group <- as.factor(sim_dat$group)

outcome_vars <- c("x1", "x2", "x3", "x4")

## Just check if function throws any errors
test_that("2 treatments, no conf_int, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
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
                                            method="direct_pseudo",
                                            conf_int=F,
                                            bootstrap=T,
                                            n_boot=5,
                                            outcome_vars=outcome_vars,
                                            type_time="bs"), NA)
})

test_that("2 treatments, no conf_int, no boot, with times", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct_pseudo",
                                            conf_int=F,
                                            bootstrap=T,
                                            n_boot=5,
                                            outcome_vars=outcome_vars,
                                            times=c(1, 2),
                                            type_time="bs"), NA)
})

sim_dat$time <- round(sim_dat$time, 1)

test_that("2 treatments, no conf_int, no boot, with times, type_time='factor'", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct_pseudo",
                                            conf_int=F,
                                            bootstrap=T,
                                            n_boot=5,
                                            outcome_vars=outcome_vars,
                                            times=c(1, 2),
                                            type_time="factor"), NA)
})

sim_dat <- adjustedCurves::sim_confounded_surv(n=300)
sim_dat$group[sim_dat$group==1] <- sample(c(1, 2),
                                          size=nrow(sim_dat[sim_dat$group==1,]),
                                          replace=T)
sim_dat$group <- as.factor(sim_dat$group)
sim_dat$time <- round(sim_dat$time, 1)

test_that("> 2 treatments, no conf_int, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
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
                                            method="direct_pseudo",
                                            conf_int=F,
                                            bootstrap=T,
                                            n_boot=5,
                                            outcome_vars=outcome_vars,
                                            type_time="bs"), NA)
})

test_that("> 2 treatments, no conf_int, no boot, with times", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct_pseudo",
                                            conf_int=F,
                                            bootstrap=F,
                                            n_boot=5,
                                            outcome_vars=outcome_vars,
                                            times=c(1, 2),
                                            type_time="bs"), NA)
})

test_that("> 2 treatments, no conf_int, no boot, with times, type_time='factor'", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct_pseudo",
                                            conf_int=F,
                                            bootstrap=F,
                                            n_boot=5,
                                            outcome_vars=outcome_vars,
                                            times=c(1, 2),
                                            type_time="factor"), NA)
})

