library(survival)

set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_surv_n_50.Rds",
                               package="adjustedCurves"))
sim_dat$group <- as.factor(sim_dat$group)

# outcome model
outcome_vars <- c("x2", "x3", "x4", "x6")

# treatment model
treat_mod <- glm(group ~ x1 + x2 + x3 + x4 + x5 + x6, data=sim_dat,
                 family="binomial")
ps_score <- stats::predict(treat_mod, newdata=sim_dat, type="response")

test_that("2 treatments, no conf_int, no boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="aiptw_pseudo",
                      conf_int=FALSE,
                      outcome_vars=outcome_vars,
                      treatment_model=treat_mod,
                      type_time="bs")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("2 treatments, with conf_int, no boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="aiptw_pseudo",
                      conf_int=TRUE,
                      outcome_vars=outcome_vars,
                      treatment_model=treat_mod,
                      type_time="bs")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("2 treatments, no conf_int, with boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="aiptw_pseudo",
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=2,
                      outcome_vars=outcome_vars,
                      treatment_model=treat_mod,
                      type_time="bs")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("2 treatments, no conf_int, no boot, with times, factor time", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="aiptw_pseudo",
                      conf_int=FALSE,
                      bootstrap=FALSE,
                      n_boot=2,
                      outcome_vars=outcome_vars,
                      treatment_model=treat_mod,
                      times=c(0.5, 0.8, 1),
                      type_time="factor")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("2 treatments, no conf_int, no boot, with times, ns time", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="aiptw_pseudo",
                      conf_int=FALSE,
                      bootstrap=FALSE,
                      n_boot=2,
                      outcome_vars=outcome_vars,
                      treatment_model=treat_mod,
                      times=c(0.5, 0.8, 1),
                      type_time="ns",
                      spline_df=2)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("2 treatments, with ps_scores supplied", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="aiptw_pseudo",
                      conf_int=FALSE,
                      bootstrap=FALSE,
                      outcome_vars=outcome_vars,
                      treatment_model=ps_score,
                      type_time="bs")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_surv_n_50.Rds",
                               package="adjustedCurves"))
sim_dat$group[sim_dat$group==1] <- sample(c(1, 2),
                                        size=nrow(sim_dat[sim_dat$group==1, ]),
                                        replace=TRUE)
sim_dat$group <- as.factor(sim_dat$group)

treat_mod <- quiet(nnet::multinom(group ~ x1 + x2, data=sim_dat))

test_that("> 2 treatments, no conf_int, no boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="aiptw_pseudo",
                      conf_int=FALSE,
                      outcome_vars=outcome_vars,
                      treatment_model=treat_mod,
                      type_time="bs")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("> 2 treatments, with conf_int, no boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="aiptw_pseudo",
                      conf_int=TRUE,
                      outcome_vars=outcome_vars,
                      treatment_model=treat_mod,
                      type_time="bs")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("> 2 treatments, no conf_int, with boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="aiptw_pseudo",
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=2,
                      outcome_vars=outcome_vars,
                      treatment_model=treat_mod,
                      type_time="bs")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("> 2 treatments, no conf_int, no boot, with times, factor time", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="aiptw_pseudo",
                      conf_int=FALSE,
                      bootstrap=FALSE,
                      n_boot=2,
                      outcome_vars=outcome_vars,
                      treatment_model=treat_mod,
                      times=c(0.5, 0.8, 1),
                      type_time="factor")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})
