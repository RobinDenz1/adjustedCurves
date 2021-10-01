library(survival)
library(riskRegression)
library(mice)

set.seed(42)

sim_dat <- adjustedCurves::sim_confounded_crisk(n=150)
sim_dat$group_num <- sim_dat$group
sim_dat$group <- as.factor(sim_dat$group)
sim_dat$x1 <- ifelse(runif(n=150) <= 0.4, sim_dat$x1, NA)

# impute dataset
imp <- mice::mice(sim_dat, m=3, method="pmm", printFlag=F)

# outcome model
outc_mod <- CSC_MI(mids=imp,
                   formula=Hist(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6)

# treatment model
treat_mod <- with(imp, glm(group ~ x1 + x2 + x3 + x4 + x5 + x6, family="binomial"))

# censoring model
cens_mod <- with(imp, coxph(Surv(time, event==0) ~ x1 + x2, x=T))

outc_vars <- c("x1", "x2", "x3", "x4", "x5", "x6")
treat_vars <- c("x1", "x2", "x3", "x4", "x5", "x6")

### direct
test_that("MI, direct, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=imp,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=F,
                                           outcome_model=outc_mod,
                                           cause=1), NA)
})

test_that("MI, direct, boot", {
  expect_error(adjustedCurves::adjustedcif(data=imp,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct",
                                           conf_int=F,
                                           bootstrap=T,
                                           n_boot=3,
                                           outcome_model=outc_mod,
                                           cause=1), NA)
})

### direct_pseudo
test_that("MI, direct_pseudo, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=imp,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct_pseudo",
                                           conf_int=F,
                                           outcome_vars=outc_vars,
                                           type_time="bs",
                                           cause=1), NA)
})

test_that("MI, direct_pseudo, boot", {
  expect_error(adjustedCurves::adjustedcif(data=imp,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="direct_pseudo",
                                           conf_int=F,
                                           bootstrap=T,
                                           n_boot=3,
                                           outcome_vars=outc_vars,
                                           type_time="bs",
                                           cause=1), NA)
})

### iptw
# use glm treatment model
test_that("MI, iptw, no boot, glm", {
  expect_error(adjustedCurves::adjustedcif(data=imp,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="iptw",
                                           conf_int=F,
                                           treatment_model=treat_mod,
                                           cause=1), NA)
})

test_that("MI, iptw, boot, glm", {
  expect_error(adjustedCurves::adjustedcif(data=imp,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="iptw",
                                           conf_int=F,
                                           bootstrap=T,
                                           n_boot=3,
                                           treatment_model=treat_mod,
                                           cause=1), NA)
})

### iptw_pseudo
# use glm model
test_that("MI, iptw_pseudo, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=imp,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="iptw_pseudo",
                                           conf_int=F,
                                           treatment_model=treat_mod,
                                           cause=1), NA)
})

test_that("MI, iptw_pseudo, boot", {
  expect_error(adjustedCurves::adjustedcif(data=imp,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="iptw_pseudo",
                                           conf_int=F,
                                           bootstrap=T,
                                           n_boot=3,
                                           treatment_model=treat_mod,
                                           cause=1), NA)
})

# use formula object treatment model
test_that("MI, iptw_pseudo, no boot, weightit", {
  expect_error(adjustedCurves::adjustedcif(data=imp,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="iptw_pseudo",
                                           conf_int=F,
                                           treatment_model=group ~ x1 + x2 + x3,
                                           cause=1), NA)
})

test_that("MI, iptw_pseudo, boot, weightit", {
  expect_error(adjustedCurves::adjustedcif(data=imp,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="iptw_pseudo",
                                           conf_int=F,
                                           bootstrap=T,
                                           n_boot=3,
                                           treatment_model=group ~ x1 + x2 + x3,
                                           cause=1), NA)
})

### matching
test_that("MI, matching, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=imp,
                                           variable="group_num",
                                           ev_time="time",
                                           event="event",
                                           method="matching",
                                           conf_int=F,
                                           treatment_model=treat_mod,
                                           cause=1), NA)
})

### aiptw
test_that("MI, aiptw, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=imp,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aiptw",
                                           conf_int=F,
                                           outcome_model=outc_mod,
                                           treatment_model=treat_mod,
                                           cause=1), NA)
})

test_that("MI, aiptw, boot", {
  expect_error(adjustedCurves::adjustedcif(data=imp,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aiptw",
                                           conf_int=F,
                                           bootstrap=T,
                                           n_boot=3,
                                           outcome_model=outc_mod,
                                           treatment_model=treat_mod,
                                           cause=1), NA)
})

### aiptw_pseudo
test_that("MI, aiptw_pseudo, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=imp,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aiptw_pseudo",
                                           conf_int=F,
                                           outcome_vars=outc_vars,
                                           treatment_model=treat_mod,
                                           type_time="bs",
                                           cause=1), NA)
})

test_that("MI, aiptw_pseudo, boot", {
  expect_error(adjustedCurves::adjustedcif(data=imp,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aiptw_pseudo",
                                           conf_int=F,
                                           bootstrap=T,
                                           n_boot=3,
                                           outcome_vars=outc_vars,
                                           treatment_model=treat_mod,
                                           type_time="bs",
                                           cause=1), NA)
})

### aalen_johansen
test_that("MI, aalen_johansen, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=imp,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aalen_johansen",
                                           conf_int=F,
                                           cause=1), NA)
})

test_that("MI, aalen_johansen, boot", {
  expect_error(adjustedCurves::adjustedcif(data=imp,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aalen_johansen",
                                           conf_int=F,
                                           bootstrap=T,
                                           n_boot=3,
                                           cause=1), NA)
})

### test_curve_equality
adjcif <- adjustedCurves::adjustedcif(data=imp,
                                      variable="group",
                                      ev_time="time",
                                      event="event",
                                      method="aalen_johansen",
                                      bootstrap=T,
                                      n_boot=3,
                                      na.action="na.omit",
                                      cause=1)

test_that("test_curve_equality, two treatments", {
  expect_error(adjustedCurves::test_curve_equality(adjcif, from=0, to=1), NA)
})

# create 3 treatments
sim_dat$group2 <- 0
sim_dat$group2[sim_dat$group==1] <- sample(c(1, 2), size=nrow(sim_dat[sim_dat$group==1,]),
                                           replace=T)
sim_dat$group2 <- ifelse(sim_dat$group2==1, "Placebo", ifelse(sim_dat$group2==2, "Chemo", "OP"))
sim_dat$group2 <- factor(sim_dat$group2)

imp <- mice::mice(sim_dat, m=3, method="pmm", printFlag=F)

# fit adjustedsurv
adjcif <- adjustedCurves::adjustedcif(data=imp,
                                      variable="group2",
                                      ev_time="time",
                                      event="event",
                                      method="aalen_johansen",
                                      bootstrap=T,
                                      n_boot=3,
                                      na.action="na.omit",
                                      cause=1)

test_that("test_curve_equality, three treatments", {
  expect_error(adjustedCurves::test_curve_equality(adjcif, from=0, to=1), NA)
})

