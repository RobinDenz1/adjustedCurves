library(survival)
library(prodlim)

set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_surv_n_20.Rds",
                               package="adjustedCurves"))
sim_dat$group <- as.factor(sim_dat$group)

sim_dat_tibble <- dplyr::tibble(sim_dat)

# adjustedsurv works with tibbles?
test_that("adjustedsurv, tibbles, no boot", {
  adj <- adjustedsurv(data=sim_dat_tibble,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_km",
                      conf_int=FALSE,
                      treatment_model=group ~ x1)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("adjustedsurv, tibbles, with boot", {
  adj <- adjustedsurv(data=sim_dat_tibble,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_km",
                      bootstrap=TRUE,
                      n_boot=2,
                      conf_int=FALSE,
                      treatment_model=group ~ x1,
                      na.action=stats::na.omit)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

sim_dat2 <- readRDS(system.file("testdata",
                                "d_sim_crisk_n_20.Rds",
                                package="adjustedCurves"))
sim_dat2$group <- as.factor(sim_dat2$group)

sim_dat_tibble <- dplyr::tibble(sim_dat2)

# adjustedcif works with tibbles?
test_that("adjustedcif, tibbles, no boot", {
  adj <- adjustedcif(data=sim_dat_tibble,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="aalen_johansen",
                     conf_int=FALSE,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("adjustedcif, tibbles, with boot", {
  adj <- adjustedcif(data=sim_dat_tibble,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="aalen_johansen",
                     bootstrap=TRUE,
                     n_boot=2,
                     conf_int=FALSE,
                     cause=1,
                     na.action=stats::na.omit)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

## multicore processing

test_that("adjustedsurv, km, with n_cores = 2", {
  adj <- quiet(adjustedsurv(data=sim_dat,
                            variable="group",
                            ev_time="time",
                            event="event",
                            method="km",
                            bootstrap=TRUE,
                            n_boot=2,
                            conf_int=FALSE,
                            n_cores=2))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("adjustedsurv, aiptw, with n_cores = 2", {
  cox_mod <- survival::coxph(Surv(time, event) ~ group + x1,
                             data=sim_dat, x=TRUE)
  ps_mod <- glm(group ~ x1, data=sim_dat, family="binomial")
  cens_mod <- survival::coxph(Surv(time, event==0) ~ x2 + x3,
                              data=sim_dat, x=TRUE)

  adj <- quiet(adjustedsurv(data=sim_dat,
                            variable="group",
                            ev_time="time",
                            event="event",
                            method="aiptw",
                            bootstrap=TRUE,
                            n_boot=2,
                            conf_int=FALSE,
                            n_cores=2,
                            outcome_model=cox_mod,
                            censoring_model=cens_mod,
                            treatment_model=ps_mod))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("adjustedcif, aalen_johansen, with n_cores = 2", {
  adj <- quiet(adjustedcif(data=sim_dat2,
                           variable="group",
                           ev_time="time",
                           event="event",
                           method="aalen_johansen",
                           bootstrap=TRUE,
                           n_boot=2,
                           conf_int=FALSE,
                           cause=1,
                           n_cores=2))
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("adjustedcif, aiptw, with n_cores = 2", {
  csc_mod <- riskRegression::CSC(Hist(time, event) ~ group + x4, data=sim_dat2)
  ps_mod <- glm(group ~ x1, data=sim_dat2, family="binomial")
  cens_mod <- survival::coxph(Surv(time, event==0) ~ x2 + x3,
                              data=sim_dat2, x=TRUE)

  adj <- quiet(suppressWarnings(adjustedcif(data=sim_dat2,
                                variable="group",
                                ev_time="time",
                                event="event",
                                method="aiptw",
                                bootstrap=TRUE,
                                n_boot=2,
                                conf_int=FALSE,
                                cause=1,
                                n_cores=2,
                                outcome_model=csc_mod,
                                censoring_model=cens_mod,
                                treatment_model=ps_mod)))
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

## warnings and errors

test_that("adjustedsurv, data wrong format", {
  expect_error(adjustedsurv(data="?",
                            variable="group",
                            ev_time="time",
                            event="event",
                            method="km",
                            bootstrap=FALSE,
                            conf_int=FALSE),
               "'data' must be either a data.frame or mids object.")
})

test_that("adjustedcif, data wrong format", {
  expect_error(adjustedcif(data="?",
                           variable="group",
                           ev_time="time",
                           event="event",
                           method="km",
                           bootstrap=FALSE,
                           conf_int=FALSE),
               "'data' must be either a data.frame or mids object.")
})

test_that("adjustedsurv, all data removed by na.action call", {
  sim_dat_err <- sim_dat
  sim_dat_err$x1 <- NA
  expect_error(adjustedsurv(data=sim_dat_err,
                            variable="group",
                            ev_time="time",
                            event="event",
                            method="iptw_km",
                            treatment_model=group ~ x1,
                            bootstrap=FALSE,
                            conf_int=FALSE,
                            na.action="na.omit"),
               "There is no non-missing data left after call to 'na.action'.")
})

test_that("adjustedcif, all data removed by na.action call", {
  sim_dat_err <- sim_dat2
  sim_dat_err$x1 <- NA
  expect_error(adjustedcif(data=sim_dat_err,
                           variable="group",
                           ev_time="time",
                           event="event",
                           method="iptw_pseudo",
                           treatment_model=group ~ x1,
                           bootstrap=FALSE,
                           conf_int=FALSE,
                           na.action="na.omit",
                           cause=1),
               "There is no non-missing data left after call to 'na.action'.")
})

test_that("adjustedsurv, NA in relevant data when using weights", {
  sim_dat_err <- sim_dat
  sim_dat_err$group[1] <- NA
  weights <- runif(n=nrow(sim_dat_err), min=1, max=3)

  expect_error(adjustedsurv(data=sim_dat_err,
                            variable="group",
                            ev_time="time",
                            event="event",
                            method="iptw_pseudo",
                            treatment_model=weights,
                            bootstrap=FALSE,
                            conf_int=FALSE,
                            na.action="na.omit"),
               paste0("Weights cannot be supplied directly to the",
                      " 'treatment_model' argument if there are missing",
                      " values in relevant columns of 'data'."))
})

test_that("adjustedcif, NA in relevant data when using weights", {
  sim_dat_err <- sim_dat2
  sim_dat_err$group[1] <- NA
  weights <- runif(n=nrow(sim_dat_err), min=1, max=3)

  expect_error(adjustedcif(data=sim_dat_err,
                           variable="group",
                           ev_time="time",
                           event="event",
                           method="iptw_pseudo",
                           treatment_model=weights,
                           bootstrap=FALSE,
                           conf_int=FALSE,
                           na.action="na.omit",
                           cause=1),
               paste0("Weights cannot be supplied directly to the",
                      " 'treatment_model' argument if there are missing",
                      " values in relevant columns of 'data'."))
})

test_that("adjustedsurv, using force_bounds and iso_reg with bootstrapping", {
  set.seed(435)
  adj_without <- adjustedsurv(data=sim_dat,
                              variable="group",
                              ev_time="time",
                              event="event",
                              method="iptw_pseudo",
                              treatment_model=group ~ x1 + x2,
                              bootstrap=TRUE,
                              n_boot=10,
                              force_bounds=FALSE,
                              iso_reg=FALSE,
                              conf_int=FALSE,
                              na.action="na.omit")
  min_without <- min(adj_without$boot_data$surv, na.rm=TRUE)
  max_without <- max(adj_without$boot_data$surv, na.rm=TRUE)

  set.seed(435)
  adj_with <- adjustedsurv(data=sim_dat,
                           variable="group",
                           ev_time="time",
                           event="event",
                           method="iptw_pseudo",
                           treatment_model=group ~ x1 + x2,
                           bootstrap=TRUE,
                           n_boot=10,
                           force_bounds=TRUE,
                           iso_reg=TRUE,
                           conf_int=FALSE,
                           na.action="na.omit")
  min_with <- min(adj_with$boot_data$surv, na.rm=TRUE)
  max_with <- max(adj_with$boot_data$surv, na.rm=TRUE)

  expect_true(min_without < 0)
  expect_true(max_without > 1)
  expect_true(min_with==0)
  expect_true(max_with==1)
})

test_that("adjustedcif, using force_bounds and iso_reg with bootstrapping", {
  set.seed(435)
  adj_without <- suppressWarnings(adjustedcif(data=sim_dat,
                              variable="group",
                              ev_time="time",
                              event="event",
                              method="iptw_pseudo",
                              treatment_model=group ~ x1 + x2,
                              bootstrap=TRUE,
                              n_boot=10,
                              force_bounds=FALSE,
                              iso_reg=FALSE,
                              conf_int=FALSE,
                              na.action="na.omit",
                             cause=1))
  min_without <- min(adj_without$boot_data$cif, na.rm=TRUE)
  max_without <- max(adj_without$boot_data$cif, na.rm=TRUE)

  set.seed(435)
  adj_with <- suppressWarnings(adjustedcif(data=sim_dat,
                           variable="group",
                           ev_time="time",
                           event="event",
                           method="iptw_pseudo",
                           treatment_model=group ~ x1 + x2,
                           bootstrap=TRUE,
                           n_boot=10,
                           force_bounds=TRUE,
                           iso_reg=TRUE,
                           conf_int=FALSE,
                           na.action="na.omit",
                           cause=1))
  min_with <- min(adj_with$boot_data$cif, na.rm=TRUE)
  max_with <- max(adj_with$boot_data$cif, na.rm=TRUE)

  expect_true(min_without < 0)
  expect_true(max_without > 1)
  expect_true(min_with==0)
  expect_true(max_with==1)
})
