library(prodlim)

set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_crisk_n_100.Rds",
                               package="adjustedCurves"))
sim_dat$group <- as.factor(sim_dat$group)

# outcome models
mod_CSC <- riskRegression::CSC(Hist(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6,
                           data=sim_dat)
mod_FGR <- riskRegression::FGR(Hist(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6,
                            data=sim_dat, cause=1)
# fit models
mod_riskRegression <- riskRegression::riskRegression(Hist(time, event) ~ x1 +
                                                       x2 + x3 + group,
                                                     data=sim_dat, cause=1)
mod_ARR <- riskRegression::ARR(Hist(time, event) ~ x1 + x2 + x3 + group,
                               data=sim_dat, cause=1)
mod_prodlim <- prodlim::prodlim(prodlim::Hist(time, event) ~ x1 + x2 + x3 +
                                  group, data=sim_dat)

## Just check if function throws any errors
test_that("CSC, 2 treatments, no conf_int, no boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=FALSE,
                     outcome_model=mod_CSC,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("CSC, 2 treatments, with conf_int, no boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=TRUE,
                     outcome_model=mod_CSC,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("CSC, 2 treatments, no conf_int, with boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=FALSE,
                     bootstrap=TRUE,
                     n_boot=2,
                     outcome_model=mod_CSC,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("CSC, 2 treatments, with conf_int, with boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=TRUE,
                     bootstrap=TRUE,
                     n_boot=2,
                     outcome_model=mod_CSC,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("CSC, 2 treatments, no conf_int, no boot, with times", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=TRUE,
                     bootstrap=TRUE,
                     n_boot=2,
                     outcome_model=mod_CSC,
                     times=c(0.5, 1, 1.2),
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("FGR, 2 treatments, no boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=FALSE,
                     outcome_model=mod_FGR,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("FGR, 2 treatments, with boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=FALSE,
                     bootstrap=TRUE,
                     n_boot=2,
                     outcome_model=mod_FGR,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("FGR, 2 treatments, no boot, with times", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=FALSE,
                     bootstrap=FALSE,
                     n_boot=2,
                     outcome_model=mod_FGR,
                     times=c(0.5, 1, 1.2),
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("FGR, 2 treatments, with boot, with times", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=FALSE,
                     bootstrap=TRUE,
                     n_boot=2,
                     outcome_model=mod_FGR,
                     times=c(0.5, 1, 1.2),
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("riskRegression, 2 treatments, no boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=FALSE,
                     outcome_model=mod_riskRegression,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("riskRegression, 2 treatments, with boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     bootstrap=TRUE,
                     n_boot=2,
                     outcome_model=mod_riskRegression,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("ARR, 2 treatments, no boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=FALSE,
                     outcome_model=mod_ARR,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("ARR, 2 treatments, with boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     bootstrap=TRUE,
                     n_boot=2,
                     outcome_model=mod_ARR,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("prodlim, 2 treatments, no boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=FALSE,
                     outcome_model=mod_prodlim,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

mod_fastcmprsk <- readRDS(system.file("testdata",
                                      "fastcmprsk_model.Rds",
                                      package="adjustedCurves"))

test_that("fastcmprsk, 2 treatments, no boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=FALSE,
                     outcome_model=mod_fastcmprsk,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

set.seed(37)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_crisk_n_200.Rds",
                               package="adjustedCurves"))
sim_dat$group[sim_dat$group==1] <- sample(c(1, 2),
                                        size=nrow(sim_dat[sim_dat$group==1, ]),
                                        replace=TRUE)
sim_dat$group <- as.factor(sim_dat$group)

# outcome models
mod_CSC <- riskRegression::CSC(Hist(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6,
                           data=sim_dat)
mod_FGR <- riskRegression::FGR(Hist(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6,
                            data=sim_dat, cause=1)
mod_riskRegression <- riskRegression::riskRegression(Hist(time, event) ~ x1 +
                                                       x2 + x3 + group,
                                                     data=sim_dat, cause=1)
mod_ARR <- riskRegression::ARR(Hist(time, event) ~ x1 + x2 + x3 + group,
                               data=sim_dat, cause=1)
mod_prodlim <- prodlim::prodlim(prodlim::Hist(time, event) ~ x1 + x2 + x3 +
                                  group, data=sim_dat)

test_that("CSC, > 2 treatments, no conf_int, no boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=FALSE,
                     outcome_model=mod_CSC,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("CSC, > 2 treatments, with conf_int, no boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=TRUE,
                     outcome_model=mod_CSC,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("CSC, > 2 treatments, no conf_int, with boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=FALSE,
                     bootstrap=TRUE,
                     n_boot=2,
                     outcome_model=mod_CSC,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("CSC, > 2 treatments, with conf_int, with boot", {
  adj <- suppressWarnings(adjustedcif(data=sim_dat,
                                      variable="group",
                                      ev_time="time",
                                      event="event",
                                      method="direct",
                                      conf_int=TRUE,
                                      bootstrap=TRUE,
                                      n_boot=2,
                                      outcome_model=mod_CSC,
                                      cause=1))
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("CSC, > 2 treatments, no conf_int, no boot, with times", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=FALSE,
                     bootstrap=FALSE,
                     n_boot=2,
                     outcome_model=mod_CSC,
                     times=c(0.5, 1, 1.2),
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("CSC, > 2 treatments, no conf_int, with boot, with times", {
  adj <- suppressWarnings(adjustedcif(data=sim_dat,
                          variable="group",
                          ev_time="time",
                          event="event",
                          method="direct",
                          conf_int=FALSE,
                          bootstrap=TRUE,
                          n_boot=2,
                          outcome_model=mod_CSC,
                          times=c(0.5, 1, 1.2),
                          cause=1))
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("FGR, > 2 treatments, no boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=FALSE,
                     outcome_model=mod_FGR,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("FGR, > 2 treatments, with boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=FALSE,
                     bootstrap=TRUE,
                     n_boot=2,
                     outcome_model=mod_FGR,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("FGR, > 2 treatments, no boot, with times", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=FALSE,
                     bootstrap=FALSE,
                     n_boot=2,
                     outcome_model=mod_FGR,
                     times=c(0.5, 1, 1.2),
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("FGR, > 2 treatments, with boot, with times", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=FALSE,
                     bootstrap=TRUE,
                     n_boot=2,
                     outcome_model=mod_FGR,
                     times=c(0.5, 1, 1.2),
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("riskRegression, > 2 treatments, no boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=FALSE,
                     outcome_model=mod_riskRegression,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("riskRegression, > 2 treatments, with boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=TRUE,
                     outcome_model=mod_riskRegression,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("ARR, > 2 treatments, no boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=FALSE,
                     outcome_model=mod_ARR,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("ARR, > 2 treatments, with boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=TRUE,
                     outcome_model=mod_ARR,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("prodlim, > 2 treatments, no boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=FALSE,
                     outcome_model=mod_prodlim,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("prodlim, > 2 treatments, with boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=TRUE,
                     outcome_model=mod_prodlim,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

## models from the mets package
if (requireNamespace("mets", quietly=TRUE)) {

  sim_dat_2 <- readRDS(system.file("testdata",
                                   "d_sim_crisk_n_100.Rds",
                                   package="adjustedCurves"))
  sim_dat_2$group <- as.factor(sim_dat_2$group)

  mod_cifreg_2 <- mets::cifregFG(mets::Event(time, event) ~ group + x1 + x2 +
                                   x3, data=sim_dat_2, cause=1,
                                 control=list(tol=1e-10, iter.max=100))
  mod_cifreg_3 <- mets::cifregFG(mets::Event(time, event) ~ x1 + x2 + x3 +
                                   group, data=sim_dat, cause=1)
}

test_that("cifreg, 2 treatments, no conf_int", {
  skip_if_not_installed("mets")
  adj <- adjustedcif(data=sim_dat_2,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=FALSE,
                     outcome_model=mod_cifreg_2,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat_2$group))
})

test_that("cifreg, 2 treatments, with conf_int", {
  skip_if_not_installed("mets")
  adj <- adjustedcif(data=sim_dat_2,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=TRUE,
                     outcome_model=mod_cifreg_2,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$se))
  expect_true(all(adj$adj$ci_lower <= adj$adj$cif))
  expect_true(all(adj$adj$ci_upper >= adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat_2$group))
})

test_that("cifreg, 2 treatments, with times", {
  skip_if_not_installed("mets")
  adj <- adjustedcif(data=sim_dat_2,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=TRUE,
                     outcome_model=mod_cifreg_2,
                     times=c(0, 0.5, 1, 1.2),
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_equal(nrow(adj$adj), 8)
  expect_equal(adj$adj$cif[adj$adj$time==0], c(0, 0))
  expect_equal(levels(adj$adj$group), levels(sim_dat_2$group))
})

test_that("cifreg, > 2 treatments, with conf_int", {
  skip_if_not_installed("mets")
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     conf_int=TRUE,
                     outcome_model=mod_cifreg_3,
                     times=c(0.5, 1, 1.2),
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$se))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("cifreg, same estimates as FGR with the same model", {
  skip_if_not_installed("mets")
  mod_FGR_2 <- riskRegression::FGR(Hist(time, event) ~ group + x1 + x2 + x3,
                                   data=sim_dat_2, cause=1)
  times <- c(0.2, 0.5, 1, 1.2)
  adj_mets <- adjustedcif(data=sim_dat_2,
                          variable="group",
                          ev_time="time",
                          event="event",
                          method="direct",
                          conf_int=FALSE,
                          outcome_model=mod_cifreg_2,
                          times=times,
                          cause=1)
  adj_fgr <- adjustedcif(data=sim_dat_2,
                         variable="group",
                         ev_time="time",
                         event="event",
                         method="direct",
                         conf_int=FALSE,
                         outcome_model=mod_FGR_2,
                         times=times,
                         cause=1)
  # crr() itself only converges to about 1e-5
  expect_equal(adj_mets$adj$cif, adj_fgr$adj$cif, tolerance=1e-4)
})

test_that("cifreg, all tied events are included at a tied event time", {
  skip_if_not_installed("mets")
  sim_dat_tied <- sim_dat_2
  sim_dat_tied$time <- ceiling(sim_dat_tied$time * 10)
  mod_tied <- mets::cifregFG(mets::Event(time, event) ~ group + x1 + x2 + x3,
                             data=sim_dat_tied, cause=1)

  # a time with several events, no observed times in (t, t + 0.5]
  n_events <- table(sim_dat_tied$time[sim_dat_tied$event==1])
  t_tied <- as.numeric(names(n_events)[n_events > 1][1])

  for (conf_int in c(FALSE, TRUE)) {
    adj <- adjustedcif(data=sim_dat_tied,
                       variable="group",
                       ev_time="time",
                       event="event",
                       method="direct",
                       conf_int=conf_int,
                       outcome_model=mod_tied,
                       times=c(t_tied, t_tied + 0.5),
                       cause=1)
    expect_equal(adj$adj$cif[adj$adj$time==t_tied],
                 adj$adj$cif[adj$adj$time==t_tied + 0.5])
  }
})

test_that("cifreg, same estimates with and without conf_int", {
  skip_if_not_installed("mets")
  times <- c(1, 0, 0.2, 0.5)
  adj_ci <- adjustedcif(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="direct",
                        conf_int=TRUE,
                        outcome_model=mod_cifreg_3,
                        times=times,
                        cause=1)
  adj_no_ci <- adjustedcif(data=sim_dat,
                           variable="group",
                           ev_time="time",
                           event="event",
                           method="direct",
                           conf_int=FALSE,
                           outcome_model=mod_cifreg_3,
                           times=times,
                           cause=1)
  expect_equal(adj_ci$adj$time, adj_no_ci$adj$time)
  expect_equal(adj_ci$adj$cif, adj_no_ci$adj$cif)
})

# NOTE: These models are supported, but would require all kinds of package
#       dependencies if tested here. The tests are run in the development
#       process but are commented out here.

#library(randomForestSRC)
#mod_rfsrc <- randomForestSRC::rfsrc(Surv(time, event) ~ x1 + x2 + x3 + group,
#                                   data=sim_dat)
#library(casebase)
#mod_fitSmoothHazard <- casebase::fitSmoothHazard(event ~ time + x1 + group,
#                                                 sim_dat, ratio=10)
#
#
#test_that("rfsrc, 2 treatments, no boot", {
#  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
#                                           variable="group",
#                                           ev_time="time",
#                                           event="event",
#                                           method="direct",
#                                           conf_int=FALSE,
#                                           outcome_model=mod_rfsrc,
#                                           cause=1), NA)
#})
#
#test_that("rfsrc, 2 treatments, with boot", {
#  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
#                                           variable="group",
#                                           ev_time="time",
#                                           event="event",
#                                           method="direct",
#                                           conf_int=TRUE,
#                                           outcome_model=mod_rfsrc,
#                                           cause=1), NA)
#})
#
## NOTE: Currently doesn't work due to bugs in casebase
#test_that("fitSmoothHazard, 2 treatments, no boot", {
#  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
#                                           variable="group",
#                                           ev_time="time",
#                                           event="event",
#                                           method="direct",
#                                           conf_int=FALSE,
#                                           outcome_model=mod_fitSmoothHazard,
#                                           cause=1), NA)
#})
#
## NOTE Currently doesn't work due to bugs in casebase
#test_that("fitSmoothHazard, 2 treatments, with boot", {
#  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
#                                           variable="group",
#                                           ev_time="time",
#                                           event="event",
#                                           method="direct",
#                                           conf_int=TRUE,
#                                           outcome_model=mod_fitSmoothHazard,
#                                           cause=1), NA)
#})
