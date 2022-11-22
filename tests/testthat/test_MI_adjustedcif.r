library(survival)

set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_crisk_n_50.Rds",
                               package="adjustedCurves"))

sim_dat_tibble <- dplyr::tibble(sim_dat)
sim_dat$group <- as.factor(sim_dat$group)
sim_dat$x1 <- ifelse(runif(n=nrow(sim_dat)) <= 0.7, sim_dat$x1, NA)

# impute dataset
imp <- suppressWarnings(mice::mice(sim_dat, m=3, method="pmm", printFlag=FALSE))

# outcome model
outc_mod <- CSC_MI(mids=imp,
                   formula=Hist(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6)

# treatment model
treat_mod <- with(imp, glm(group ~ x1 + x2 + x3 + x4 + x5 + x6,
                           family="binomial"))

# censoring model
cens_mod <- with(imp, coxph(Surv(time, event==0) ~ x1 + x2, x=TRUE))

outc_vars <- c("x1", "x2", "x3", "x4", "x5", "x6")
treat_vars <- c("x1", "x2", "x3", "x4", "x5", "x6")

### direct
test_that("MI, direct, no boot", {
  adj <- adjustedcif(data=imp,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=FALSE,
                     outcome_model=outc_mod,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

test_that("MI, direct, boot", {
  suppressWarnings({
  adj <- adjustedcif(data=imp,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=FALSE,
                     bootstrap=TRUE,
                     n_boot=2,
                     outcome_model=outc_mod,
                     cause=1)})
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

test_that("MI, direct, conf_int", {
  adj <- adjustedcif(data=imp,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=TRUE,
                     bootstrap=FALSE,
                     n_boot=2,
                     outcome_model=outc_mod,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

### direct_pseudo
test_that("MI, direct_pseudo, no boot", {
  adj <- adjustedcif(data=imp,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct_pseudo",
                     conf_int=FALSE,
                     outcome_vars=outc_vars,
                     type_time="bs",
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

test_that("MI, direct_pseudo, boot", {
  adj <- adjustedcif(data=imp,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct_pseudo",
                     conf_int=FALSE,
                     bootstrap=TRUE,
                     n_boot=2,
                     outcome_vars=outc_vars,
                     type_time="bs",
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

### iptw
# use glm treatment model
test_that("MI, iptw, no boot, glm", {
  adj <- adjustedcif(data=imp,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="iptw",
                     conf_int=FALSE,
                     treatment_model=treat_mod,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

test_that("MI, iptw, boot, glm", {
  suppressWarnings({
  adj <- adjustedcif(data=imp,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="iptw",
                     conf_int=FALSE,
                     bootstrap=TRUE,
                     n_boot=2,
                     treatment_model=treat_mod,
                     cause=1)})
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

### iptw_pseudo
# use glm model
test_that("MI, iptw_pseudo, no boot", {
  adj <- adjustedcif(data=imp,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="iptw_pseudo",
                     conf_int=FALSE,
                     treatment_model=treat_mod,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

test_that("MI, iptw_pseudo, boot", {
  adj <- suppressWarnings({adjustedcif(data=imp,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="iptw_pseudo",
                     conf_int=FALSE,
                     bootstrap=TRUE,
                     n_boot=2,
                     treatment_model=treat_mod,
                     cause=1)})
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

# use formula object treatment model
test_that("MI, iptw_pseudo, no boot, weightit", {
  adj <- adjustedcif(data=imp,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="iptw_pseudo",
                     conf_int=FALSE,
                     treatment_model=group ~ x1 + x2 + x3,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

test_that("MI, iptw_pseudo, boot, weightit", {
  adj <- adjustedcif(data=imp,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="iptw_pseudo",
                     conf_int=FALSE,
                     bootstrap=TRUE,
                     n_boot=2,
                     treatment_model=group ~ x1 + x2 + x3,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

### matching
test_that("MI, matching, no boot", {
  adj <- adjustedcif(data=imp,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="matching",
                     conf_int=FALSE,
                     treatment_model=treat_mod,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

### aiptw
test_that("MI, aiptw, no boot", {
  adj <- adjustedcif(data=imp,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="aiptw",
                     conf_int=FALSE,
                     outcome_model=outc_mod,
                     treatment_model=treat_mod,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

test_that("MI, aiptw, boot", {
  adj <- suppressWarnings({adjustedcif(data=imp,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="aiptw",
                     conf_int=FALSE,
                     bootstrap=TRUE,
                     n_boot=2,
                     outcome_model=outc_mod,
                     treatment_model=treat_mod,
                     cause=1)})
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

### aiptw_pseudo
test_that("MI, aiptw_pseudo, no boot", {
  adj <- adjustedcif(data=imp,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="aiptw_pseudo",
                     conf_int=FALSE,
                     outcome_vars=outc_vars,
                     treatment_model=treat_mod,
                     type_time="bs",
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

test_that("MI, aiptw_pseudo, boot", {
  adj <- adjustedcif(data=imp,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="aiptw_pseudo",
                     conf_int=FALSE,
                     bootstrap=TRUE,
                     n_boot=2,
                     outcome_vars=outc_vars,
                     treatment_model=treat_mod,
                     type_time="bs",
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

### aalen_johansen
test_that("MI, aalen_johansen, no boot", {
  adj <- adjustedcif(data=imp,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="aalen_johansen",
                     conf_int=F,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

test_that("MI, aalen_johansen, boot", {
  adj <- adjustedcif(data=imp,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="aalen_johansen",
                     conf_int=F,
                     bootstrap=T,
                     n_boot=2,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

### adjusted_curve_test
adjcif <- adjustedcif(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="aalen_johansen",
                      bootstrap=T,
                      n_boot=2,
                      na.action="na.omit",
                      cause=1)

test_that("adjusted_curve_test, two treatments", {
  adj_test <- adjusted_curve_test(adjcif, from=0, to=1)
  expect_equal(round(adj_test$observed_diff_integral, 4), 0.0129)
  expect_equal(round(adj_test$integral_se, 3), 0.053)
  expect_equal(adj_test$mids_p_values, c(1, 1, 1))
  expect_equal(adj_test$n_boot, 2)
})

test_that("adjusted_curve_test, two treatments, linear", {
  adj_test <- adjusted_curve_test(adjcif, from=0, to=1, interpolation="linear")
  expect_equal(round(adj_test$observed_diff_integral, 4), 0.0135)
  expect_equal(round(adj_test$integral_se, 3), 0.053)
  expect_equal(adj_test$mids_p_values, c(1, 1, 1))
  expect_equal(adj_test$n_boot, 2)
})

# create 3 treatments
sim_dat$group2 <- 0
sim_dat$group2[sim_dat$group==1] <-
  sample(c(1, 2), size=nrow(sim_dat[sim_dat$group==1, ]), replace=T)
sim_dat$group2 <- ifelse(sim_dat$group2==1, "Placebo",
                         ifelse(sim_dat$group2==2, "Chemo", "OP"))
sim_dat$group2 <- factor(sim_dat$group2)

imp <- suppressWarnings(mice::mice(sim_dat, m=3, method="pmm", printFlag=F))

# fit adjustedsurv
adjcif <- adjustedcif(data=imp,
                      variable="group2",
                      ev_time="time",
                      event="event",
                      method="aalen_johansen",
                      bootstrap=T,
                      n_boot=2,
                      na.action="na.omit",
                      cause=1)

test_that("adjusted_curve_test, three treatments", {
  adj_test <- adjusted_curve_test(adjcif, from=0, to=1)
  expect_equal(round(adj_test$`Chemo vs. Placebo`$observed_diff_integral, 1),
               -0.1)
  expect_true(is.na(adj_test$`Chemo vs. Placebo`$integral_se))
  expect_equal(round(adj_test$`Chemo vs. Placebo`$p_value, 4), 0)
  expect_equal(adj_test$`Chemo vs. Placebo`$mids_p_values, c(0, 0, 0))
  expect_equal(adj_test$`Chemo vs. Placebo`$n_boot, 1)
})

test_that("adjusted_curve_test, three treatments linear", {
  adj_test <- adjusted_curve_test(adjcif, from=0, to=1, interpolation="linear")
  expect_equal(round(adj_test$`Chemo vs. Placebo`$observed_diff_integral, 1),
               -0.1)
  expect_true(is.na(adj_test$`Chemo vs. Placebo`$integral_se))
  expect_equal(round(adj_test$`Chemo vs. Placebo`$p_value, 4), 0)
  expect_equal(adj_test$`Chemo vs. Placebo`$mids_p_values, c(0, 0, 0))
  expect_equal(adj_test$`Chemo vs. Placebo`$n_boot, 1)
})
