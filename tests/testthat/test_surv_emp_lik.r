library(survival)

set.seed(42)

sim_dat <- adjustedCurves::sim_confounded_surv(n=100)
sim_dat <- within(sim_dat, {
  x1 <- ifelse(x1==1, 1, -1)
  x2 <- ifelse(x2==1, 1, -1)
  x3 <- ifelse(x3==1, 1, -1)
})

# outcome model
treatment_vars <- c("x1", "x2", "x3", "x4", "x5")

## Just check if function throws any errors
test_that("2 treatments, no conf_int, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="emp_lik",
                                            conf_int=F,
                                            treatment_vars=treatment_vars), NA)
})

test_that("2 treatments, no conf_int, with boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="emp_lik",
                                            conf_int=F,
                                            bootstrap=T,
                                            n_boot=2,
                                            treatment_vars=treatment_vars), NA)
})

test_that("2 treatments, no conf_int, no boot, with times", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="emp_lik",
                                            conf_int=F,
                                            bootstrap=F,
                                            treatment_vars=treatment_vars,
                                            times=c(0.5, 0.8, 1)), NA)
})
