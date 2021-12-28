library(survival)

set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_surv_n_50.Rds",
                               package="adjustedCurves"))
sim_dat$group <- as.factor(sim_dat$group)

# outcome model
outc_mod <- survival::coxph(Surv(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6,
                            data=sim_dat, x=TRUE)

# treatment model
treat_mod <- glm(group ~ x1 + x2 + x3 + x4 + x5 + x6, data=sim_dat,
                 family="binomial")

## Just check if function throws any errors
test_that("2 treatments, no conf_int, no boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="aiptw",
                      conf_int=FALSE,
                      outcome_model=outc_mod,
                      treatment_model=treat_mod)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("2 treatments, with conf_int, no boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="aiptw",
                      conf_int=TRUE,
                      outcome_model=outc_mod,
                      treatment_model=treat_mod)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("2 treatments, no conf_int, with boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="aiptw",
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=2,
                      outcome_model=outc_mod,
                      treatment_model=treat_mod)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("2 treatments, no conf_int, no boot, with times", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="aiptw",
                      conf_int=FALSE,
                      bootstrap=FALSE,
                      n_boot=2,
                      outcome_model=outc_mod,
                      treatment_model=treat_mod,
                      times=c(0.8, 0.9))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

cens_mod <- survival::coxph(Surv(time, event==0) ~ x2, data=sim_dat, x=TRUE)

test_that("2 treatments, no conf_int, no boot, with times, with cens_mod", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="aiptw",
                      conf_int=FALSE,
                      bootstrap=FALSE,
                      n_boot=2,
                      outcome_model=outc_mod,
                      treatment_model=treat_mod,
                      censoring_model=cens_mod,
                      times=c(0.8, 0.9))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})
