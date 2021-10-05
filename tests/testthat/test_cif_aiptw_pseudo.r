library(adjustedCurves)
library(survival)
library(riskRegression)

set.seed(42)

sim_dat <- adjustedCurves::sim_confounded_surv(n=50)
sim_dat$event[sim_dat$event==1] <- sample(c(1, 2), size=sum(sim_dat$event),
                                          replace=T)
sim_dat$group <- as.factor(sim_dat$group)

# outcome model
outcome_vars <- c("x2", "x3", "x4", "x6")

# treatment model
treat_mod <- glm(group ~ x1 + x2 + x3 + x4 + x5 + x6, data=sim_dat,
                 family="binomial")

## Just check if function throws any errors
test_that("2 treatments, no conf_int, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aiptw_pseudo",
                                           conf_int=F,
                                           outcome_vars=outcome_vars,
                                           treatment_model=treat_mod,
                                           type_time="bs",
                                           cause=1), NA)
})

test_that("2 treatments, with conf_int, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aiptw_pseudo",
                                           conf_int=T,
                                           outcome_vars=outcome_vars,
                                           treatment_model=treat_mod,
                                           type_time="bs",
                                           cause=1), NA)
})

test_that("2 treatments, no conf_int, with boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aiptw_pseudo",
                                           conf_int=F,
                                           bootstrap=T,
                                           n_boot=2,
                                           outcome_vars=outcome_vars,
                                           treatment_model=treat_mod,
                                           type_time="bs",
                                           cause=1), NA)
})

test_that("2 treatments, no conf_int, no boot, with times, factor time", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aiptw_pseudo",
                                           conf_int=F,
                                           bootstrap=F,
                                           n_boot=2,
                                           outcome_vars=outcome_vars,
                                           treatment_model=treat_mod,
                                           times=c(0.5, 0.8, 1),
                                           type_time="factor",
                                           cause=1), NA)
})


sim_dat <- adjustedCurves::sim_confounded_surv(n=50)
sim_dat$group[sim_dat$group==1] <- sample(c(1, 2),
                                          size=nrow(sim_dat[sim_dat$group==1,]),
                                          replace=T)
sim_dat$group <- as.factor(sim_dat$group)
sim_dat$event[sim_dat$event==1] <- sample(c(1, 2), size=sum(sim_dat$event),
                                          replace=T)

treat_mod <- nnet::multinom(group ~ x1 + x2, data=sim_dat)

test_that("> 2 treatments, no conf_int, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aiptw_pseudo",
                                           conf_int=F,
                                           outcome_vars=outcome_vars,
                                           treatment_model=treat_mod,
                                           type_time="bs",
                                           cause=1), NA)
})

test_that("> 2 treatments, with conf_int, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aiptw_pseudo",
                                           conf_int=T,
                                           outcome_vars=outcome_vars,
                                           treatment_model=treat_mod,
                                           type_time="bs",
                                           cause=1), NA)
})

test_that("> 2 treatments, no conf_int, with boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aiptw_pseudo",
                                           conf_int=F,
                                           bootstrap=T,
                                           n_boot=2,
                                           outcome_vars=outcome_vars,
                                           treatment_model=treat_mod,
                                           type_time="bs",
                                           cause=1), NA)
})

test_that("> 2 treatments, no conf_int, no boot, with times, factor time", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aiptw_pseudo",
                                           conf_int=F,
                                           bootstrap=F,
                                           n_boot=2,
                                           outcome_vars=outcome_vars,
                                           treatment_model=treat_mod,
                                           times=c(0.5, 0.8, 1),
                                           type_time="factor",
                                           cause=1), NA)
})
