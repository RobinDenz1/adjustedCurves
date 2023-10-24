
set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_surv_n_50.Rds",
                               package="adjustedCurves"))
sim_dat$group <- ifelse(sim_dat$group==0, "Control", "Treatment")
sim_dat$group <- as.factor(sim_dat$group)

test_that("one model each", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="tmle",
                      conf_int=TRUE,
                      treatment_model=c("SL.glm"),
                      outcome_model=list(Surv(time, status) ~ .),
                      times=c(0.8, 1))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("multiple models each", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="tmle",
                      conf_int=TRUE,
                      treatment_model=c("SL.glm", "SL.mean"),
                      outcome_model=list(Surv(time, status) ~ x1 + x2,
                                         Surv(time, status) ~ x3 + x5),
                      times=c(0.8, 1))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("no conf_int", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="tmle",
                      conf_int=FALSE,
                      treatment_model=c("SL.glm"),
                      outcome_model=list(Surv(time, status) ~ .),
                      times=c(0.8, 1))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
  expect_true(length(colnames(adj$adjsurv)) == 3)
})

test_that("changing some arguments", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="tmle",
                      conf_int=TRUE,
                      treatment_model=c("SL.glm"),
                      outcome_model=list(Surv(time, status) ~ .),
                      times=c(0.8, 1),
                      cv_args=list(V=2),
                      max_update_iter=150,
                      one_step_eps=0.12,
                      min_nuisance=0.1,
                      return_models=FALSE)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})
