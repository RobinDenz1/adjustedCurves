library(survival)
library(riskRegression)
library(mice)

set.seed(42)

sim_dat <- adjustedCurves::sim_confounded_surv(n=100)
sim_dat$group_num <- sim_dat$group
sim_dat$group <- as.factor(sim_dat$group)
sim_dat$x1 <- ifelse(runif(n=100) <= 0.4, sim_dat$x1, NA)

# impute dataset
imp <- mice::mice(sim_dat, m=3, method="pmm", printFlag=FALSE)

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
  expect_error(adjustedCurves::adjustedsurv(data=imp,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            outcome_model=outc_mod), NA)
})

test_that("MI, direct, boot", {
  expect_error(adjustedCurves::adjustedsurv(data=imp,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            bootstrap=TRUE,
                                            n_boot=3,
                                            outcome_model=outc_mod), NA)
})

### direct_pseudo
test_that("MI, direct_pseudo, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=imp,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct_pseudo",
                                            conf_int=FALSE,
                                            outcome_vars=outc_vars,
                                            type_time="bs"), NA)
})

test_that("MI, direct_pseudo, boot", {
  expect_error(adjustedCurves::adjustedsurv(data=imp,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct_pseudo",
                                            conf_int=FALSE,
                                            bootstrap=TRUE,
                                            n_boot=3,
                                            outcome_vars=outc_vars,
                                            type_time="bs"), NA)
})

### iptw_km
# use glm treatment model
test_that("MI, iptw_km, no boot, glm", {
  expect_error(adjustedCurves::adjustedsurv(data=imp,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="iptw_km",
                                            conf_int=FALSE,
                                            treatment_model=treat_mod), NA)
})

test_that("MI, iptw_km, boot, glm", {
  expect_error(adjustedCurves::adjustedsurv(data=imp,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="iptw_km",
                                            conf_int=FALSE,
                                            bootstrap=TRUE,
                                            n_boot=3,
                                            treatment_model=treat_mod), NA)
})

# use formula object treatment model
test_that("MI, iptw_km, no boot, weightit", {
  expect_error(adjustedCurves::adjustedsurv(data=imp,
                                          variable="group",
                                          ev_time="time",
                                          event="event",
                                          method="iptw_km",
                                          conf_int=FALSE,
                                          treatment_model=group ~ x1 + x2 + x3),
               NA)
})

test_that("MI, iptw_km, boot, weightit", {
  expect_error(adjustedCurves::adjustedsurv(data=imp,
                                          variable="group",
                                          ev_time="time",
                                          event="event",
                                          method="iptw_km",
                                          conf_int=FALSE,
                                          bootstrap=TRUE,
                                          n_boot=3,
                                          treatment_model=group ~ x1 + x2 + x3),
               NA)
})


### iptw_cox
# use glm model
test_that("MI, iptw_cox, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=imp,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="iptw_cox",
                                            conf_int=FALSE,
                                            treatment_model=treat_mod), NA)
})

test_that("MI, iptw_cox, boot", {
  expect_error(adjustedCurves::adjustedsurv(data=imp,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="iptw_cox",
                                            conf_int=FALSE,
                                            bootstrap=TRUE,
                                            n_boot=3,
                                            treatment_model=treat_mod), NA)
})

# use formula object treatment model
test_that("MI, iptw_cox, no boot, weightit", {
  expect_error(adjustedCurves::adjustedsurv(data=imp,
                                          variable="group",
                                          ev_time="time",
                                          event="event",
                                          method="iptw_cox",
                                          conf_int=FALSE,
                                          treatment_model=group ~ x1 + x2 + x3),
               NA)
})

test_that("MI, iptw_cox, boot, weightit", {
  expect_error(adjustedCurves::adjustedsurv(data=imp,
                                          variable="group",
                                          ev_time="time",
                                          event="event",
                                          method="iptw_cox",
                                          conf_int=FALSE,
                                          bootstrap=TRUE,
                                          n_boot=3,
                                          treatment_model=group ~ x1 + x2 + x3),
               NA)
})

### iptw_pseudo
# use glm model
test_that("MI, iptw_pseudo, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=imp,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="iptw_pseudo",
                                            conf_int=FALSE,
                                            treatment_model=treat_mod), NA)
})

test_that("MI, iptw_pseudo, boot", {
  expect_error(adjustedCurves::adjustedsurv(data=imp,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="iptw_pseudo",
                                            conf_int=FALSE,
                                            bootstrap=TRUE,
                                            n_boot=3,
                                            treatment_model=treat_mod), NA)
})

# use formula object treatment model
test_that("MI, iptw_pseudo, no boot, weightit", {
  expect_error(adjustedCurves::adjustedsurv(data=imp,
                                          variable="group",
                                          ev_time="time",
                                          event="event",
                                          method="iptw_pseudo",
                                          conf_int=FALSE,
                                          treatment_model=group ~ x1 + x2 + x3),
               NA)
})

test_that("MI, iptw_pseudo, boot, weightit", {
  expect_error(adjustedCurves::adjustedsurv(data=imp,
                                          variable="group",
                                          ev_time="time",
                                          event="event",
                                          method="iptw_pseudo",
                                          conf_int=FALSE,
                                          bootstrap=TRUE,
                                          n_boot=3,
                                          treatment_model=group ~ x1 + x2 + x3),
               NA)
})

### matching
test_that("MI, matching, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=imp,
                                            variable="group_num",
                                            ev_time="time",
                                            event="event",
                                            method="matching",
                                            conf_int=FALSE,
                                            treatment_model=treat_mod), NA)
})

### aiptw
test_that("MI, aiptw, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=imp,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="aiptw",
                                            conf_int=FALSE,
                                            outcome_model=outc_mod,
                                            treatment_model=treat_mod), NA)
})

test_that("MI, aiptw, boot", {
  expect_error(adjustedCurves::adjustedsurv(data=imp,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="aiptw",
                                            conf_int=FALSE,
                                            bootstrap=TRUE,
                                            n_boot=3,
                                            outcome_model=outc_mod,
                                            treatment_model=treat_mod), NA)
})

### aiptw_pseudo
test_that("MI, aiptw_pseudo, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=imp,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="aiptw_pseudo",
                                            conf_int=FALSE,
                                            outcome_vars=outc_vars,
                                            treatment_model=treat_mod,
                                            type_time="bs"), NA)
})

test_that("MI, aiptw_pseudo, boot", {
  expect_error(adjustedCurves::adjustedsurv(data=imp,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="aiptw_pseudo",
                                            conf_int=FALSE,
                                            bootstrap=TRUE,
                                            n_boot=3,
                                            outcome_vars=outc_vars,
                                            treatment_model=treat_mod,
                                            type_time="bs"), NA)
})

### km
test_that("MI, km, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=imp,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="km",
                                            conf_int=FALSE), NA)
})

test_that("MI, km, boot", {
  expect_error(adjustedCurves::adjustedsurv(data=imp,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="km",
                                            conf_int=FALSE,
                                            bootstrap=TRUE,
                                            n_boot=3), NA)
})

### emp_lik
sim_dat$x1 <- ifelse(sim_dat$x1==0, -1, 1)
sim_dat$x2 <- ifelse(sim_dat$x2==0, -1, 1)
sim_dat$x3 <- ifelse(sim_dat$x3==0, -1, 1)

imp <- mice::mice(sim_dat, m=3, method="pmm", printFlag=FALSE)

test_that("MI, emp_lik, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=imp,
                                            variable="group_num",
                                            ev_time="time",
                                            event="event",
                                            method="emp_lik",
                                            conf_int=FALSE,
                                            treatment_vars=treat_vars), NA)
})

test_that("MI, emp_lik, boot", {
  expect_error(adjustedCurves::adjustedsurv(data=imp,
                                            variable="group_num",
                                            ev_time="time",
                                            event="event",
                                            method="emp_lik",
                                            conf_int=FALSE,
                                            bootstrap=TRUE,
                                            n_boot=3,
                                            treatment_vars=treat_vars), NA)
})

### adjusted_median_survival
adjsurv <- adjustedCurves::adjustedsurv(data=imp,
                                        variable="group",
                                        ev_time="time",
                                        event="event",
                                        method="km",
                                        bootstrap=TRUE,
                                        n_boot=3,
                                        na.action="na.omit")

test_that("adjusted_median_survival, no boot", {
  expect_error(adjustedCurves::adjusted_median_survival(adjsurv,
                                                        use_boot=FALSE),
               NA)
})

test_that("adjusted_median_survival, boot", {
  expect_error(adjustedCurves::adjusted_median_survival(adjsurv,
                                                        use_boot=TRUE),
               NA)
})

### adjusted_rmst
test_that("adjusted_rmst, no boot", {
  expect_error(adjustedCurves::adjusted_rmst(adjsurv, from=0, to=1,
                                             use_boot=FALSE),
               NA)
})

test_that("adjusted_rmst, boot", {
  expect_error(adjustedCurves::adjusted_rmst(adjsurv, from=0, to=1,
                                             use_boot=TRUE),
               NA)
})

### adjusted_curve_diff
test_that("adjusted_curve_diff, two treatments", {
  expect_error(adjustedCurves::adjusted_curve_diff(adjsurv, from=0, to=1), NA)
})

# create 3 treatments
sim_dat$group2 <- 0
sim_dat$group2[sim_dat$group==1] <-
  sample(c(1, 2), size=nrow(sim_dat[sim_dat$group==1, ]), replace=TRUE)
sim_dat$group2 <- ifelse(sim_dat$group2==1, "Placebo",
                         ifelse(sim_dat$group2==2, "Chemo", "OP"))
sim_dat$group2 <- factor(sim_dat$group2)

imp <- mice::mice(sim_dat, m=3, method="pmm", printFlag=FALSE)

# fit adjustedsurv
adjsurv <- adjustedCurves::adjustedsurv(data=imp,
                                        variable="group2",
                                        ev_time="time",
                                        event="event",
                                        method="km",
                                        bootstrap=TRUE,
                                        n_boot=3,
                                        na.action="na.omit")

test_that("adjusted_curve_diff, three treatments", {
  expect_error(adjustedCurves::adjusted_curve_diff(adjsurv, from=0, to=1), NA)
})
