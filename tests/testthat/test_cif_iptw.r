library(nnet)
library(survival)
library(riskRegression)

set.seed(42)

sim_dat <- adjustedCurves::sim_confounded_surv(n=100)
sim_dat$event[sim_dat$event==1] <- sample(c(1, 2), size=sum(sim_dat$event),
                                          replace=T)
sim_dat$group <- as.factor(sim_dat$group)

mod <- glm(group ~ x1 + x2 + x3 + x4 + x5 + x6, data=sim_dat,
           family="binomial")

## Just check if function throws any errors
test_that("2 treatments, no conf_int, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="iptw",
                                           conf_int=F,
                                           treatment_model=mod,
                                           cause=1), NA)
})

test_that("2 treatments, with conf_int, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="iptw",
                                           conf_int=T,
                                           treatment_model=mod,
                                           cause=1), NA)
})

test_that("2 treatments, no conf_int, with boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="iptw",
                                           conf_int=F,
                                           bootstrap=T,
                                           n_boot=2,
                                           treatment_model=mod,
                                           cause=1), NA)
})

test_that("2 treatments, with conf_int, with boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="iptw",
                                           conf_int=T,
                                           bootstrap=T,
                                           n_boot=2,
                                           treatment_model=mod,
                                           cause=1), NA)
})

test_that("2 treatments, no conf_int, with WeightIt", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="iptw",
                                           conf_int=F,
                                           bootstrap=F,
                                           treatment_model=group ~ x1 + x2,
                                           weight_method="ps",
                                           cause=1), NA)
})

sim_dat <- adjustedCurves::sim_confounded_surv(n=150)
sim_dat$event[sim_dat$event==1] <- sample(c(1, 2), size=sum(sim_dat$event),
                                          replace=T)
sim_dat$group[sim_dat$group==1] <- sample(c(1, 2),
                                          size=nrow(sim_dat[sim_dat$group==1,]),
                                          replace=T)
sim_dat$group <- as.factor(sim_dat$group)

mod <- nnet::multinom(group ~ x1 + x2 + x3 + x4 + x5 + x6, data=sim_dat)

test_that("> 2 treatments, no conf_int, no boot, no ...", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="iptw",
                                           conf_int=F,
                                           treatment_model=mod,
                                           cause=1), NA)
})

test_that("> 2 treatments, with conf_int, no boot, no ...", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="iptw",
                                           conf_int=T,
                                           treatment_model=mod,
                                           cause=1), NA)
})

test_that("> 2 treatments, no conf_int, with boot, no ...", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="iptw",
                                           conf_int=F,
                                           bootstrap=T,
                                           n_boot=2,
                                           treatment_model=mod,
                                           cause=1), NA)
})

test_that("> 2 treatments, with conf_int, with boot, no ...", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="iptw",
                                           conf_int=T,
                                           bootstrap=T,
                                           n_boot=2,
                                           treatment_model=mod,
                                           cause=1), NA)
})

# DONT RUN: would require package dependency on "mlogit"
#test_that("> 2 treatments, no conf_int, with WeightIt", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="iptw_cox",
#                                            conf_int=F,
#                                            bootstrap=F,
#                                            treatment_model=group ~ x1 + x2,
#                                            weight_method="ps"), NA)
#})