library(survival)

set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_surv_n_150.Rds",
                               package="adjustedCurves"))
sim_dat$group_num <- sim_dat$group
sim_dat$group <- as.factor(sim_dat$group)
sim_dat$x1 <- ifelse(runif(n=nrow(sim_dat)) <= 0.7, sim_dat$x1, NA)

# impute dataset
imp <- suppressWarnings(mice::mice(sim_dat, m=3, method="pmm", printFlag=FALSE))

# outcome model
outc_mod <- with(imp, coxph(Surv(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6 +
                              group,
                            x=TRUE))

# treatment model
treat_mod <- with(imp, glm(group ~ x1 + x2 + x3 + x4 + x5 + x6,
                           family="binomial"))

# censoring model
cens_mod <- with(imp, coxph(Surv(time, event==0) ~ x1 + x2, x=TRUE))

outc_vars <- c("x1", "x2", "x3", "x4", "x5", "x6")
treat_vars <- c("x1", "x2", "x3", "x4", "x5", "x6")

### direct
test_that("MI, direct, no boot", {
  adj <- adjustedsurv(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="direct",
                      conf_int=FALSE,
                      outcome_model=outc_mod)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("MI, direct, boot", {
  adj <- adjustedsurv(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="direct",
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=2,
                      outcome_model=outc_mod)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

### direct_pseudo
test_that("MI, direct_pseudo, no boot", {
  adj <- adjustedsurv(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="direct_pseudo",
                      conf_int=FALSE,
                      outcome_vars=outc_vars,
                      type_time="bs")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("MI, direct_pseudo, boot", {
  adj <- adjustedsurv(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="direct_pseudo",
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=2,
                      outcome_vars=outc_vars,
                      type_time="bs")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

### iptw_km
# use glm treatment model
test_that("MI, iptw_km, no boot, glm", {
  adj <- adjustedsurv(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_km",
                      conf_int=FALSE,
                      treatment_model=treat_mod)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("MI, iptw_km, boot, glm", {
  adj <- adjustedsurv(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_km",
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=2,
                      treatment_model=treat_mod)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

# use formula object treatment model
test_that("MI, iptw_km, no boot, weightit", {
  adj <- adjustedsurv(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_km",
                      conf_int=FALSE,
                      treatment_model=group ~ x1 + x2 + x3)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("MI, iptw_km, boot, weightit", {
  adj <- adjustedsurv(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_km",
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=2,
                      treatment_model=group ~ x1 + x2 + x3)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})


### iptw_cox
# use glm model
test_that("MI, iptw_cox, no boot", {
  adj <- adjustedsurv(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_cox",
                      conf_int=FALSE,
                      treatment_model=treat_mod)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("MI, iptw_cox, boot", {
  adj <- adjustedsurv(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_cox",
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=2,
                      treatment_model=treat_mod)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

# use formula object treatment model
test_that("MI, iptw_cox, no boot, weightit", {
  adj <- adjustedsurv(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_cox",
                      conf_int=FALSE,
                      treatment_model=group ~ x1 + x2 + x3)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("MI, iptw_cox, boot, weightit", {
  adj <- adjustedsurv(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_cox",
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=2,
                      treatment_model=group ~ x1 + x2 + x3)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

### iptw_pseudo
# use glm model
test_that("MI, iptw_pseudo, no boot", {
  adj <- adjustedsurv(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_pseudo",
                      conf_int=FALSE,
                      treatment_model=treat_mod)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("MI, iptw_pseudo, boot", {
  adj <- adjustedsurv(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_pseudo",
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=2,
                      treatment_model=treat_mod)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

# use formula object treatment model
test_that("MI, iptw_pseudo, no boot, weightit", {
  adj <- adjustedsurv(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_pseudo",
                      conf_int=FALSE,
                      treatment_model=group ~ x1 + x2 + x3)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("MI, iptw_pseudo, boot, weightit", {
  adj <- adjustedsurv(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_pseudo",
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=2,
                      treatment_model=group ~ x1 + x2 + x3)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("MI, iptw_pseudo, using conf_int", {
  adj <- adjustedsurv(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_pseudo",
                      conf_int=TRUE,
                      bootstrap=FALSE,
                      n_boot=2,
                      treatment_model=group ~ x1 + x2 + x3)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

### matching
test_that("MI, matching, no boot", {
  adj <- adjustedsurv(data=imp,
                      variable="group_num",
                      ev_time="time",
                      event="event",
                      method="matching",
                      conf_int=FALSE,
                      treatment_model=treat_mod)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

### aiptw
test_that("MI, aiptw, no boot", {
  adj <- adjustedsurv(data=imp,
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

test_that("MI, aiptw, boot", {
  adj <- adjustedsurv(data=imp,
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

### aiptw_pseudo
test_that("MI, aiptw_pseudo, no boot", {
  adj <- adjustedsurv(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="aiptw_pseudo",
                      conf_int=FALSE,
                      outcome_vars=outc_vars,
                      treatment_model=treat_mod,
                      type_time="bs")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("MI, aiptw_pseudo, boot", {
  adj <- adjustedsurv(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="aiptw_pseudo",
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=2,
                      outcome_vars=outc_vars,
                      treatment_model=treat_mod,
                      type_time="bs")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

### km
test_that("MI, km, no boot", {
  adj <- adjustedsurv(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="km",
                      conf_int=FALSE)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("MI, km, boot", {
  adj <- adjustedsurv(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="km",
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=3)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

### emp_lik
sim_dat$x1 <- ifelse(sim_dat$x1==0, -1, 1)
sim_dat$x2 <- ifelse(sim_dat$x2==0, -1, 1)
sim_dat$x3 <- ifelse(sim_dat$x3==0, -1, 1)

imp <- suppressWarnings(mice::mice(sim_dat, m=3, method="pmm", printFlag=FALSE))

test_that("MI, emp_lik, no boot", {
  adj <- adjustedsurv(data=imp,
                      variable="group_num",
                      ev_time="time",
                      event="event",
                      method="emp_lik",
                      conf_int=FALSE,
                      treatment_vars=treat_vars)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("MI, emp_lik, boot", {
  adj <- adjustedsurv(data=imp,
                      variable="group_num",
                      ev_time="time",
                      event="event",
                      method="emp_lik",
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=3,
                      treatment_vars=treat_vars)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

### adjusted_median_survival
adjsurv <- adjustedsurv(data=imp,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="km",
                        bootstrap=TRUE,
                        n_boot=3,
                        na.action="na.omit")

test_that("adjusted_median_survival, 2 treatments, no boot", {
  adj_med <- adjusted_median_survival(adjsurv, verbose=FALSE)
  expect_equal(round(adj_med$median_surv, 4), c(0.4798, 0.5900))
})

### adjusted_rmst
test_that("adjusted_rmst, 2 treatments, no boot", {
  adj_rmst <- adjusted_rmst(adjsurv, from=0, to=1, use_boot=FALSE)
  expect_equal(as.vector(round(adj_rmst$auc, 4)), c(0.5284, 0.6288))
})

test_that("adjusted_rmst, 2 treatments, boot", {
  adj_rmst <- adjusted_rmst(adjsurv, from=0, to=1, use_boot=TRUE)
  expect_equal(as.vector(round(adj_rmst$auc, 4)), c(0.5284, 0.6288))
  expect_equal(as.vector(round(adj_rmst$auc_se, 4)), c(0.0357, 0.0447))
  expect_equal(as.vector(adj_rmst$n_boot), c(3, 3))
})

### adjusted_curve_diff
test_that("adjusted_curve_diff, 2 treatments", {
  adj_test <- adjusted_curve_diff(adjsurv, from=0, to=1)
  expect_equal(round(adj_test$observed_diff_integral, 4), -0.1004)
  expect_equal(round(adj_test$integral_se, 4), 0.0376)
  expect_equal(round(adj_test$p_value, 4), 0)
  expect_equal(adj_test$n_boot, 3)
})

### adjusted_curve_diff
test_that("adjusted_curve_diff, 2 treatments, linear", {
  adj_test <- adjusted_curve_diff(adjsurv, from=0, to=1, interpolation="linear")
  expect_equal(round(adj_test$observed_diff_integral, 4), -0.1009)
  expect_equal(round(adj_test$integral_se, 4), 0.0374)
  expect_equal(round(adj_test$p_value, 4), 0)
  expect_equal(adj_test$n_boot, 3)
})

# create 3 treatments
sim_dat$group2 <- 0
sim_dat$group2[sim_dat$group==1] <-
  sample(c(1, 2), size=nrow(sim_dat[sim_dat$group==1, ]), replace=TRUE)
sim_dat$group2 <- ifelse(sim_dat$group2==1, "Placebo",
                         ifelse(sim_dat$group2==2, "Chemo", "OP"))
sim_dat$group2 <- factor(sim_dat$group2)

imp <- suppressWarnings(mice::mice(sim_dat, m=3, method="pmm", printFlag=FALSE))

# fit adjustedsurv
adjsurv <- adjustedsurv(data=imp,
                        variable="group2",
                        ev_time="time",
                        event="event",
                        method="km",
                        bootstrap=TRUE,
                        n_boot=3,
                        na.action="na.omit")

test_that("adjusted_median_survival, 3 treatments, no boot", {
  adj_med <- adjusted_median_survival(adjsurv, verbose=FALSE)
  expect_equal(round(adj_med$median_surv, 4), c(0.7122, 0.4798, 0.5398))
})

### adjusted_rmst
test_that("adjusted_rmst, 3 treatments, no boot", {
  adj_rmst <- adjusted_rmst(adjsurv, from=0, to=1, use_boot=FALSE)
  expect_equal(as.vector(round(adj_rmst$auc, 4)), c(0.6588, 0.5284, 0.5924))
})

test_that("adjusted_rmst, 3 treatments, with boot", {
  adj_rmst <- adjusted_rmst(adjsurv, from=0, to=1, use_boot=TRUE)
  expect_equal(as.vector(round(adj_rmst$auc, 4)), c(0.6588, 0.5284, 0.5924))
  expect_equal(as.vector(round(adj_rmst$auc_se, 4)), c(0.0346, 0.0215, 0.0328))
  expect_equal(as.vector(adj_rmst$n_boot), c(3, 3, 3))
})

test_that("adjusted_curve_diff, 3 treatments", {
  adj_test <- adjusted_curve_diff(adjsurv, from=0, to=1, conf_level=0.95)
  expect_equal(round(adj_test$`Chemo vs. OP`$observed_diff_integral, 4),
               0.1304)
  expect_equal(round(adj_test$`Chemo vs. OP`$integral_se, 4), 0.0501)
  expect_equal(round(adj_test$`Chemo vs. OP`$p_value, 4), 0)
  expect_equal(round(adj_test$`Chemo vs. OP`$mids_p_values, 4),
               c(0, 0, 0))
  expect_equal(adj_test$`Chemo vs. OP`$n_boot, 3)
  expect_snapshot_output(print(adj_test))
})
