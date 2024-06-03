
set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_surv_n_50.Rds",
                               package="adjustedCurves"))
sim_dat$group <- as.factor(sim_dat$group)

outcome_vars <- c("x1", "x2", "x3", "x4")

## Geese
test_that("2 treatments, no conf_int, no boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="direct_pseudo",
                      conf_int=FALSE,
                      outcome_vars=outcome_vars,
                      type_time="bs")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("2 treatments, no conf_int, with boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="direct_pseudo",
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=2,
                      outcome_vars=outcome_vars,
                      type_time="bs")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("2 treatments, no conf_int, no boot, with times", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="direct_pseudo",
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=2,
                      outcome_vars=outcome_vars,
                      times=c(0.3, 0.8),
                      type_time="factor")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that(
  "2 treatments, no conf_int, no boot, with times, type_time='factor'", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="direct_pseudo",
                      conf_int=FALSE,
                      bootstrap=FALSE,
                      n_boot=5,
                      outcome_vars=outcome_vars,
                      times=c(0.3, 0.8),
                      type_time="factor")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that(
  "2 treatments, no conf_int, no boot, with times, type_time='ns'", {
    adj <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="direct_pseudo",
                        conf_int=FALSE,
                        bootstrap=FALSE,
                        n_boot=5,
                        outcome_vars=outcome_vars,
                        times=c(0.3, 0.8, 0.81),
                        type_time="ns",
                        spline_df=2)
    expect_s3_class(adj, "adjustedsurv")
    expect_true(is.numeric(adj$adj$surv))
    expect_equal(levels(adj$adj$group), levels(sim_dat$group))
  })

test_that("> 2 treatments, no conf_int, no boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="direct_pseudo",
                      conf_int=FALSE,
                      outcome_vars=outcome_vars,
                      type_time="bs")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("> 2 treatments, no conf_int, with boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="direct_pseudo",
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=2,
                      outcome_vars=outcome_vars,
                      type_time="bs")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("> 2 treatments, no conf_int, no boot, with times", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="direct_pseudo",
                      conf_int=FALSE,
                      bootstrap=FALSE,
                      n_boot=2,
                      outcome_vars=outcome_vars,
                      times=c(0.3, 0.8, 1),
                      type_time="factor")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("> 2 treatments, no conf_int, no boot, with times, type_time='factor'"
          , {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="direct_pseudo",
                      conf_int=FALSE,
                      bootstrap=FALSE,
                      n_boot=2,
                      outcome_vars=outcome_vars,
                      times=c(0.3, 0.8, 1),
                      type_time="factor")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})
