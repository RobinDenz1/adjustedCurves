
set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_crisk_n_50.Rds",
                               package="adjustedCurves"))
sim_dat$time <- sim_dat$time
sim_dat$group <- ifelse(sim_dat$group==0, "Control", "Treatment")
sim_dat$group <- as.factor(sim_dat$group)

test_that("one model each", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="tmle",
                     conf_int=TRUE,
                     treatment_model=c("SL.glm"),
                     outcome_model=list(Surv(time, status==1) ~ .),
                     cause=1,
                     times=c(0.8, 1))
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("multiple models each", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="tmle",
                     conf_int=TRUE,
                     treatment_model=c("SL.glm", "SL.mean"),
                     outcome_model=list(Surv(time, status==1) ~ x1 + x2,
                                        Surv(time, status==1) ~ x3 + x5),
                     cause=1,
                     times=c(0.8, 1))
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("no conf_int", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="tmle",
                     conf_int=FALSE,
                     treatment_model=c("SL.glm"),
                     outcome_model=list(Surv(time, status==1) ~ .),
                     cause=1,
                     times=c(0.8, 1))
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
  expect_true(length(colnames(adj$adj)) == 3)
})

test_that("changing some arguments", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="tmle",
                     conf_int=TRUE,
                     treatment_model=c("SL.glm"),
                     outcome_model=list(Surv(time, status==1) ~ x5),
                     cause=1,
                     times=c(0.8, 1),
                     cv_args=list(V=2),
                     max_update_iter=150,
                     one_step_eps=0.12,
                     min_nuisance=0.1,
                     return_models=FALSE)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})
