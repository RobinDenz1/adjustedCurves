
set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_surv_n_100.Rds",
                               package="adjustedCurves"))
sim_dat <- within(sim_dat, {
  x1 <- ifelse(x1==1, 1, -1)
  x2 <- ifelse(x2==1, 1, -1)
  x3 <- ifelse(x3==1, 1, -1)
  group <- factor(group)
})

# outcome model
treatment_vars <- c("x1", "x2", "x3", "x4", "x5")

## Just check if function throws any errors
test_that("2 treatments, no conf_int, no boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="emp_lik",
                      conf_int=FALSE,
                      treatment_vars=treatment_vars)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("2 treatments, no conf_int, with boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="emp_lik",
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=2,
                      treatment_vars=treatment_vars)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("2 treatments, no conf_int, no boot, with times", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="emp_lik",
                      conf_int=FALSE,
                      bootstrap=FALSE,
                      treatment_vars=treatment_vars,
                      times=c(0.5, 0.8, 1))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("instant convergence of algorithm", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="emp_lik",
                      conf_int=FALSE,
                      bootstrap=FALSE,
                      treatment_vars=treatment_vars,
                      times=1,
                      newton_tol=1)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

# NOTE: This throws a warning because the same warning shows up twice
test_that("no convergence of algorithm", {
  warns <- capture_warnings(adjustedsurv(data=sim_dat,
                              variable="group",
                              ev_time="time",
                              event="event",
                              method="emp_lik",
                              conf_int=FALSE,
                              bootstrap=FALSE,
                              treatment_vars=treatment_vars,
                              times=1,
                              newton_tol=0.0000000000001,
                              max_iter=1))
  expect_true(all(warns==paste0("Algorithm did not converge. Increasing ",
                                "the max_iter value or decreasing the ",
                                "newton_tol value might help.")))
})

test_that("using standardize", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="emp_lik",
                      conf_int=FALSE,
                      bootstrap=FALSE,
                      treatment_vars=treatment_vars,
                      times=1,
                      standardize=TRUE)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("using the second moment", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="emp_lik",
                      conf_int=FALSE,
                      bootstrap=FALSE,
                      treatment_vars=treatment_vars,
                      times=1,
                      moment="second")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

sim_dat$event[sim_dat$group==1] <- 0

test_that("no events in one group", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="emp_lik",
                      conf_int=FALSE,
                      bootstrap=FALSE,
                      treatment_vars=treatment_vars,
                      times=1)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})
