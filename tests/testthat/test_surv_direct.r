#library(adjustedCurves)
library(survival)

set.seed(42)

sim_dat <- adjustedCurves::sim_confounded_surv(n=300)
sim_dat$group <- as.factor(sim_dat$group)

# outcome model
mod <- survival::coxph(Surv(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6,
                       data=sim_dat, x=T)

## Just check if function throws any errors
test_that("2 treatments, no conf_int, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=F,
                                            outcome_model=mod), NA)
})

test_that("2 treatments, with conf_int, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=T,
                                            outcome_model=mod), NA)
})

test_that("2 treatments, no conf_int, with boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=F,
                                            bootstrap=T,
                                            n_boot=10,
                                            outcome_model=mod), NA)
})

test_that("2 treatments, with conf_int, with boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=T,
                                            bootstrap=T,
                                            n_boot=10,
                                            outcome_model=mod), NA)
})

test_that("2 treatments, no conf_int, no boot, with times", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=T,
                                            bootstrap=T,
                                            n_boot=10,
                                            outcome_model=mod,
                                            times=c(1, 2)), NA)
})

sim_dat <- adjustedCurves::sim_confounded_surv(n=300)
sim_dat$group[sim_dat$group==1] <- sample(c(1, 2),
                                          size=nrow(sim_dat[sim_dat$group==1,]),
                                          replace=T)
sim_dat$group <- as.factor(sim_dat$group)

# outcome model
mod <- survival::coxph(Surv(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6,
                       data=sim_dat, x=T)


test_that("> 2 treatments, no conf_int, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=F,
                                            outcome_model=mod), NA)
})

test_that("> 2 treatments, with conf_int, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=T,
                                            outcome_model=mod), NA)
})

test_that("> 2 treatments, no conf_int, with boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=F,
                                            bootstrap=T,
                                            n_boot=10,
                                            outcome_model=mod), NA)
})

test_that("> 2 treatments, with conf_int, with boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=T,
                                            bootstrap=T,
                                            n_boot=10,
                                            outcome_model=mod), NA)
})

test_that("> 2 treatments, no conf_int, no boot, with times", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=F,
                                            bootstrap=F,
                                            n_boot=10,
                                            outcome_model=mod,
                                            times=c(1, 2)), NA)
})



