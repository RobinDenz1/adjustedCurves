library(survival)

set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_surv_n_50.Rds",
                               package="adjustedCurves"))
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
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
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
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
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
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
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
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
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
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
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
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
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
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
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
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
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
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
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
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
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
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
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
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
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
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
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
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
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
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
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
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
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
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

### matching
test_that("MI, matching, no boot", {
  adj <- adjustedsurv(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="matching",
                      conf_int=FALSE,
                      treatment_model=treat_mod)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
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
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("MI, aiptw, boot", {
  # NOTE: there is a warning related to the glm fit in one bootstrap sample,
  #       which can safely be ignored for this test
  adj <- suppressWarnings(adjustedsurv(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="aiptw",
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=2,
                      outcome_model=outc_mod,
                      treatment_model=treat_mod))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
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
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
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
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
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
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
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
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

### emp_lik
sim_dat$x1 <- ifelse(sim_dat$x1==0, -1, 1)
sim_dat$x2 <- ifelse(sim_dat$x2==0, -1, 1)
sim_dat$x3 <- ifelse(sim_dat$x3==0, -1, 1)

imp <- suppressWarnings(mice::mice(sim_dat, m=3, method="pmm", printFlag=FALSE))

test_that("MI, emp_lik, no boot", {
  adj <- adjustedsurv(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="emp_lik",
                      conf_int=FALSE,
                      treatment_vars=treat_vars)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("MI, emp_lik, boot", {
  adj <- suppressWarnings({adjustedsurv(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="emp_lik",
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=3,
                      treatment_vars=treat_vars)})
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

set.seed(3)

### adjusted_surv_quantile
adjsurv <- adjustedsurv(data=imp,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="km",
                        bootstrap=TRUE,
                        n_boot=3,
                        na.action="na.omit")

test_that("adjusted_surv_quantile, 2 treatments, no boot", {
  adj_med <- adjusted_surv_quantile(adjsurv)
  expect_equal(round(adj_med$q_surv, 4), c(0.4785, 0.6252))
})

test_that("adjusted_surv_quantile, 2 treatments, no boot, difference", {
  adj_med <- adjusted_surv_quantile(adjsurv, contrast="diff")
  expect_equal(round(adj_med$diff, 4), -0.1468)
})

test_that("adjusted_surv_quantile, 2 treatments, no boot, difference, p", {
  adj_med <- adjusted_surv_quantile(adjsurv, contrast="diff", p=c(0.4, 0.5))
  expect_equal(round(adj_med$diff, 4), c(-0.1350, -0.1468))
})

test_that("adjusted_surv_quantile, 2 treatments, conf_int, difference", {
  adj_med <- adjusted_surv_quantile(adjsurv, contrast="diff", conf_int=TRUE)
  expect_equal(round(adj_med$diff, 4), -0.1468)
  expect_equal(round(adj_med$se, 4), 0.1271)
  expect_equal(round(adj_med$ci_lower, 4), -0.3958)
  expect_equal(round(adj_med$p_value, 4), 0.2482)
})

test_that("adjusted_surv_quantile, 2 treatments, no boot, ratio", {
  adj_med <- adjusted_surv_quantile(adjsurv, contrast="ratio")
  expect_equal(round(adj_med$ratio, 4), 0.7653)
})

### adjusted_rmst
test_that("adjusted_rmst, 2 treatments, no boot", {
  adj_rmst <- adjusted_rmst(adjsurv, from=0, to=1, conf_int=FALSE)
  expect_equal(round(adj_rmst$rmst, 4), c(0.5115, 0.6566))
})

test_that("adjusted_rmst, 2 treatments, boot", {
  adj_rmst <- adjusted_rmst(adjsurv, from=0, to=1, conf_int=TRUE)
  expect_equal(round(adj_rmst$rmst, 4), c(0.5115, 0.6566))
  expect_equal(round(adj_rmst$se, 4), c(0.0405, 0.0391))
  expect_equal(adj_rmst$n_boot, c(3, 3))
})

### adjusted_curve_test
test_that("adjusted_curve_test, 2 treatments", {
  adj_test <- adjusted_curve_test(adjsurv, from=0, to=1)
  expect_equal(round(adj_test$observed_diff_integral, 4), -0.1451)
  expect_equal(round(adj_test$integral_se, 4), 0.0585)
  expect_equal(round(adj_test$p_value, 4), 0)
  expect_equal(adj_test$n_boot, 3)
})

### adjusted_curve_test
test_that("adjusted_curve_test, 2 treatments, linear", {
  adj_test <- adjusted_curve_test(adjsurv, from=0, to=1, interpolation="linear")
  expect_equal(round(adj_test$observed_diff_integral, 4), -0.1504)
  expect_equal(round(adj_test$integral_se, 4), 0.0651)
  expect_equal(round(adj_test$p_value, 4), 0)
  expect_equal(adj_test$n_boot, 3)
})

set.seed(4)

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

test_that("adjusted_surv_quantile, 3 treatments, no boot", {
  adj_med <- adjusted_surv_quantile(adjsurv)
  expect_equal(round(adj_med$q_surv, 4), c(0.7408, 0.4785, 0.6084))
})

### adjusted_rmst
test_that("adjusted_rmst, 3 treatments, no boot", {
  adj_rmst <- adjusted_rmst(adjsurv, from=0, to=1, conf_int=FALSE)
  expect_equal(round(adj_rmst$rmst, 4), c(0.7463, 0.5115, 0.5788))
})

test_that("adjusted_rmst, 3 treatments, with boot", {
  adj_rmst <- adjusted_rmst(adjsurv, from=0, to=0.6, conf_int=TRUE)
  expect_equal(round(adj_rmst$rmst, 4), c(0.5645, 0.4305, 0.5008))
  expect_equal(round(adj_rmst$se, 4), c(0.0148, 0.0182, 0.0149))
  expect_equal(adj_rmst$n_boot, c(3, 3, 3))
})

test_that("adjusted_curve_test, 3 treatments", {
  adj_test <- adjusted_curve_test(adjsurv, from=0, to=0.6, conf_level=0.95)
  expect_equal(round(adj_test$`Chemo vs. OP`$observed_diff_integral, 4),
               0.134)
  expect_equal(round(adj_test$`Chemo vs. OP`$integral_se, 4), 0.0253)
  expect_equal(round(adj_test$`Chemo vs. OP`$p_value, 4), 0)
  expect_equal(round(adj_test$`Chemo vs. OP`$mids_p_values, 4),
               c(0, 0, 0))
  expect_equal(adj_test$`Chemo vs. OP`$n_boot, 3)
})
