library(survival)

set.seed(42)

sim_dat <- adjustedCurves::sim_confounded_surv(n=50)
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

test_that("2 treatments, no conf_int, no boot, with times, type_time='factor'",
          {
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
                                            times=c(0.3, 0.8, 1),
                                            type_time="factor"), NA)
})

test_that("> 2 treatments, no conf_int, no boot, with times, type_time='factor'"
          , {
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
                                            times=c(0.3, 0.8, 1),
                                            type_time="factor"), NA)
})
