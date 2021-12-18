library(survival)

set.seed(42)

sim_dat <- adjustedCurves::sim_confounded_surv(n=50)
sim_dat$group <- as.factor(sim_dat$group)

# outcome model
mod <- survival::coxph(Surv(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6 + group,
                       data=sim_dat, x=TRUE)

## Just check if function throws any errors
test_that("2 treatments, no conf_int, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            outcome_model=mod), NA)
})

test_that("2 treatments, with conf_int, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=TRUE,
                                            outcome_model=mod), NA)
})

test_that("2 treatments, no conf_int, with boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            bootstrap=TRUE,
                                            n_boot=2,
                                            outcome_model=mod), NA)
})

test_that("2 treatments, with conf_int, with boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=TRUE,
                                            bootstrap=TRUE,
                                            n_boot=2,
                                            outcome_model=mod), NA)
})

test_that("2 treatments, no conf_int, no boot, with times", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=TRUE,
                                            bootstrap=TRUE,
                                            n_boot=2,
                                            outcome_model=mod,
                                            times=c(0.8, 0.9)), NA)
})

sim_dat <- adjustedCurves::sim_confounded_surv(n=90)
sim_dat$group[sim_dat$group==1] <- sample(c(1, 2),
                                        size=nrow(sim_dat[sim_dat$group==1, ]),
                                        replace=TRUE)
sim_dat$group <- as.factor(sim_dat$group)

# outcome model
mod <- survival::coxph(Surv(time, event) ~ x1 + x2 + x3 + x4 + x5 + group,
                       data=sim_dat, x=T)


test_that("> 2 treatments, no conf_int, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            outcome_model=mod), NA)
})

test_that("> 2 treatments, with conf_int, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=TRUE,
                                            outcome_model=mod), NA)
})

test_that("> 2 treatments, no conf_int, with boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            bootstrap=TRUE,
                                            n_boot=2,
                                            outcome_model=mod), NA)
})

test_that("> 2 treatments, with conf_int, with boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=TRUE,
                                            bootstrap=TRUE,
                                            n_boot=2,
                                            outcome_model=mod), NA)
})

test_that("> 2 treatments, no conf_int, no boot, with times", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            bootstrap=FALSE,
                                            n_boot=2,
                                            outcome_model=mod,
                                            times=c(0.8, 0.9)), NA)
})

####################### Models other than coxph ################################

library(riskRegression)
library(prodlim)
library(pec)

set.seed(42)

sim_dat <- adjustedCurves::sim_confounded_surv(n=110)
sim_dat$group <- as.factor(sim_dat$group)

# fit some models
mod_riskRegression <- riskRegression::riskRegression(Hist(time, event) ~
                                                       group + x1,
                                     data=sim_dat, cause=1)

mod_ARR <- riskRegression::ARR(Hist(time, event) ~ group + x1 + x6,
                               data=sim_dat, cause=1)

mod_selectCox <- pec::selectCox(survival::Surv(time, event) ~ group + x1,
                                data=sim_dat)

mod_pecRpart <- pec::pecRpart(survival::Surv(time, event) ~ group + x1,
                               data=sim_dat)

mod_prodlim <- prodlim::prodlim(Hist(time, event) ~ group + x1,
                                data=sim_dat)

mod_glm <- stats::glm(time ~ group + x1, data=sim_dat, family="gaussian")

# run tests for each model
test_that("riskRegression, 2 treatments, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            outcome_model=mod_riskRegression),
               NA)
})

test_that("riskRegression, 2 treatments, with boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            outcome_model=mod_riskRegression,
                                            bootstrap=TRUE,
                                            n_boot=2),
               NA)
})

test_that("ARR, 2 treatments, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            outcome_model=mod_ARR),
               NA)
})

test_that("ARR, 2 treatments, with boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            outcome_model=mod_ARR,
                                            bootstrap=TRUE,
                                            n_boot=2),
               NA)
})

test_that("selectCox, 2 treatments, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            outcome_model=mod_selectCox),
               NA)
})

test_that("selectCox, 2 treatments, with boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            outcome_model=mod_selectCox,
                                            bootstrap=TRUE,
                                            n_boot=2),
               NA)
})

test_that("pecRpart, 2 treatments, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            outcome_model=mod_pecRpart),
               NA)
})

test_that("prodlim, 2 treatments, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            outcome_model=mod_prodlim),
               NA)
})

test_that("prodlim, 2 treatments, with boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            outcome_model=mod_prodlim,
                                            bootstrap=TRUE,
                                            n_boot=2),
               NA)
})

sim_dat$event <- 1

test_that("glm, 2 treatments, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            outcome_model=mod_glm),
               NA)
})

test_that("glm, 2 treatments, with boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            outcome_model=mod_glm,
                                            bootstrap=TRUE,
                                            n_boot=2),
               NA)
})

### using > 2 treatments
set.seed(42)

sim_dat <- adjustedCurves::sim_confounded_surv(n=110)
sim_dat$group[sim_dat$group==1] <- sample(c(1, 2),
                                        size=nrow(sim_dat[sim_dat$group==1, ]),
                                        replace=TRUE)
sim_dat$group <- as.factor(sim_dat$group)

# fit some models
mod_riskRegression <- riskRegression::riskRegression(Hist(time, event) ~
                                                       group + x1,
                                                     data=sim_dat, cause=1)

mod_ARR <- riskRegression::ARR(Hist(time, event) ~ group + x1 + x6,
                               data=sim_dat, cause=1)

mod_selectCox <- pec::selectCox(Surv(time, event) ~ group + x1,
                                data=sim_dat)

mod_pecRpart <- pec::pecRpart(Surv(time, event) ~ group + x1,
                              data=sim_dat)

mod_prodlim <- prodlim::prodlim(Hist(time, event) ~ group + x1,
                                data=sim_dat)

mod_glm <- stats::glm(time ~ group + x1, data=sim_dat, family="gaussian")

# run tests
test_that("riskRegression, > 2 treatments, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            outcome_model=mod_riskRegression),
               NA)
})

test_that("riskRegression, > 2 treatments, with boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            outcome_model=mod_riskRegression,
                                            bootstrap=TRUE,
                                            n_boot=2),
               NA)
})

test_that("ARR, 2 treatments, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            outcome_model=mod_ARR),
               NA)
})

test_that("ARR, 2 treatments, with boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            outcome_model=mod_ARR,
                                            bootstrap=TRUE,
                                            n_boot=2),
               NA)
})

test_that("selectCox, > 2 treatments, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            outcome_model=mod_selectCox),
               NA)
})

test_that("selectCox, > 2 treatments, with boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            outcome_model=mod_selectCox,
                                            bootstrap=TRUE,
                                            n_boot=2),
               NA)
})

test_that("pecRpart, > 2 treatments, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            outcome_model=mod_pecRpart),
               NA)
})

test_that("prodlim, > 2 treatments, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            outcome_model=mod_prodlim),
               NA)
})

test_that("prodlim, > 2 treatments, with boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            outcome_model=mod_prodlim,
                                            bootstrap=TRUE,
                                            n_boot=2),
               NA)
})

sim_dat$event <- 1

test_that("glm, > 2 treatments, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            outcome_model=mod_glm),
               NA)
})

test_that("glm, > 2 treatments, with boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="direct",
                                            conf_int=FALSE,
                                            outcome_model=mod_glm,
                                            bootstrap=TRUE,
                                            n_boot=2),
               NA)
})

################################################################################
# NOTE: These models are supported, but would require all kinds of package
#       dependencies if tested here. The tests are run in the development
#       process but are commented out here.
#
#mod_pecCforest <- pec::pecCforest(Surv(time, event) ~ group + x1,
#                                  data=sim_dat,
#                                  control=party::cforest_unbiased(mtry=2))
#library(timereg)
#mod_aalen <- timereg::aalen(Surv(time, event) ~ group + x1,
#                            data=sim_dat, max.time=7, n.sim=100)
#mod_cox.aalen <- timereg::cox.aalen(Surv(time, event) ~ prop(group) + prop(x1),
#                                    data=sim_dat)
#library(rms)
#mod_psm <- rms::psm(Surv(time, event) ~ group + x1,
#                    data=sim_dat)
#mod_ols <- rms::ols(time ~ group + x1, data=sim_dat, x=TRUE, y=TRUE)
#library(flexsurv)
#mod_flexsurvreg <- flexsurv::flexsurvreg(Surv(time, event) ~ group + x1,
#                                         data=sim_dat, dist="gengamma")
#library(randomForest)
#mod_randomForest <- randomForest::randomForest(time ~ group + x1, data=sim_dat)
#library(ranger)
#mod_ranger <- ranger::ranger(Surv(time, event) ~ group + x1 + x6, data=sim_dat)
#library(randomForestSRC)
#mod_rfsrc <- randomForestSRC::rfsrc(Surv(time, event) ~ group + x1 + x6,
#                                    data=sim_dat)
#library(penalized)
#mod_penalizedS3 <- riskRegression::penalizedS3(Surv(time, event) ~ group + x1 +
#                                               x6, data=sim_dat,
#                                             lambda1=1, maxiter=1)
#library(gbm)
#mod_gbm <- gbm::gbm(Surv(time, event) ~ group + x1 + x6, data=sim_dat,
#                    distribution="coxph")
#library(casebase)
#mod_fitSmoothHazard <- casebase::fitSmoothHazard(event ~ time + x1 + group,
#                                                 sim_dat, ratio=10)
#
#library(mexhaz)
#mod_mexhaz <- mexhaz::mexhaz(Surv(time, event) ~ group + x1 + x6, data=sim_dat)
#
#test_that("pecCforest, 2 treatments, no boot", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="direct",
#                                            conf_int=FALSE,
#                                            outcome_model=mod_pecCforest),
#               NA)
#})
#
#test_that("pecCforest, 2 treatments, with boot", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="direct",
#                                            conf_int=FALSE,
#                                            outcome_model=mod_pecCforest,
#                                            bootstrap=TRUE,
#                                            n_boot=2),
#               NA)
#})
#
#test_that("aalen, 2 treatments, no boot", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="direct",
#                                            conf_int=FALSE,
#                                            outcome_model=mod_aalen),
#               NA)
#})
#
#test_that("aalen, 2 treatments, with boot", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="direct",
#                                            conf_int=FALSE,
#                                            outcome_model=mod_aalen,
#                                            bootstrap=TRUE,
#                                            n_boot=2),
#               NA)
#})
#
## NOTE: cox.aalen currently doesn't work due to bugs in predictRisk and
##       predictSurvProb
#test_that("cox.aalen, 2 treatments, no boot", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="direct",
#                                            conf_int=FALSE,
#                                            outcome_model=mod_cox.aalen),
#               NA)
#})
#
## NOTE: cox.aalen currently doesn't work due to bugs in predictRisk and
##       predictSurvProb
#test_that("cox.aalen, 2 treatments, with boot", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="direct",
#                                            conf_int=FALSE,
#                                            outcome_model=mod_cox.aalen,
#                                            bootstrap=TRUE,
#                                            n_boot=2),
#               NA)
#})
#
#test_that("psm, 2 treatments, no boot", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="direct",
#                                            conf_int=FALSE,
#                                            outcome_model=mod_psm),
#               NA)
#})
#
#test_that("psm, 2 treatments, with boot", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="direct",
#                                            conf_int=FALSE,
#                                            outcome_model=mod_psm,
#                                            bootstrap=TRUE,
#                                            n_boot=2),
#               NA)
#})
#
#test_that("flexsurvreg, 2 treatments, no boot", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="direct",
#                                            conf_int=FALSE,
#                                            outcome_model=mod_flexsurvreg),
#               NA)
#})
#
#test_that("flexsurvreg, 2 treatments, with boot", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="direct",
#                                            conf_int=FALSE,
#                                            outcome_model=mod_flexsurvreg,
#                                            bootstrap=TRUE,
#                                            n_boot=2),
#               NA)
#})
#
#test_that("randomForest, 2 treatments, no boot", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="direct",
#                                            conf_int=FALSE,
#                                            outcome_model=mod_randomForest),
#               NA)
#})
#
#test_that("randomForest, 2 treatments, with boot", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="direct",
#                                            conf_int=FALSE,
#                                            outcome_model=mod_randomForest,
#                                            bootstrap=TRUE,
#                                            n_boot=2),
#               NA)
#})
#
#test_that("ols, 2 treatments, no boot", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="direct",
#                                            conf_int=FALSE,
#                                            outcome_model=mod_ols),
#               NA)
#})
#
#test_that("ols, 2 treatments, with boot", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="direct",
#                                            conf_int=FALSE,
#                                            outcome_model=mod_ols,
#                                            bootstrap=TRUE,
#                                            n_boot=2),
#               NA)
#})
#
#test_that("ranger, 2 treatments, no boot", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="direct",
#                                            conf_int=FALSE,
#                                            outcome_model=mod_ranger),
#               NA)
#})
#
#test_that("ranger, 2 treatments, with boot", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="direct",
#                                            conf_int=FALSE,
#                                            outcome_model=mod_ranger,
#                                            bootstrap=TRUE,
#                                            n_boot=2),
#               NA)
#})
#
#test_that("rfsrc, 2 treatments, no boot", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="direct",
#                                            conf_int=FALSE,
#                                            outcome_model=mod_rfsrc),
#               NA)
#})
#
#test_that("rfsrc, 2 treatments, with boot", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="direct",
#                                            conf_int=FALSE,
#                                            outcome_model=mod_rfsrc,
#                                            bootstrap=TRUE,
#                                            n_boot=2),
#               NA)
#})
#
#test_that("penalizedS3, 2 treatments, no boot", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="direct",
#                                            conf_int=FALSE,
#                                            outcome_model=mod_penalizedS3),
#               NA)
#})
#
#test_that("penalizedS3, 2 treatments, with boot", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="direct",
#                                            conf_int=FALSE,
#                                            outcome_model=mod_penalizedS3,
#                                            bootstrap=TRUE,
#                                            n_boot=2),
#               NA)
#})
#
#test_that("gbm, 2 treatments, no boot", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="direct",
#                                            conf_int=FALSE,
#                                            outcome_model=mod_gbm),
#               NA)
#})
#
#test_that("gbm, 2 treatments, with boot", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="direct",
#                                            conf_int=FALSE,
#                                            outcome_model=mod_gbm,
#                                            bootstrap=TRUE,
#                                            n_boot=2),
#               NA)
#})
#
## NOTE: Currently doesn't work due to bugs in casebase
#test_that("fitSmoothHazard, 2 treatments, no boot", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="direct",
#                                            conf_int=FALSE,
#                                            outcome_model=mod_fitSmoothHazard),
#               NA)
#})
#
## NOTE: Currently doesn't work due to bugs in casebase
#test_that("fitSmoothHazard, 2 treatments, with boot", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="direct",
#                                            conf_int=FALSE,
#                                            outcome_model=mod_fitSmoothHazard,
#                                            bootstrap=TRUE,
#                                            n_boot=2),
#               NA)
#})
#
#test_that("mexhaz, 2 treatments, no boot", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="direct",
#                                            conf_int=FALSE,
#                                            outcome_model=mod_mexhaz),
#               NA)
#})
#
#test_that("mexhaz, 2 treatments, with boot", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="direct",
#                                            conf_int=FALSE,
#                                            outcome_model=mod_mexhaz,
#                                            bootstrap=TRUE,
#                                            n_boot=2),
#               NA)
#})
