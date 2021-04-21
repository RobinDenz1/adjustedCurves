library(survival)
library(riskRegression)

set.seed(42)

sim_dat <- adjustedCurves::sim_confounded_surv(n=100)
sim_dat$event[sim_dat$event==1] <- sample(c(1, 2), size=sum(sim_dat$event),
                                          replace=T)
sim_dat$group <- as.factor(sim_dat$group)

# outcome model
outc_mod <- riskRegression::CSC(Hist(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6,
                                data=sim_dat)

# treatment model
treat_mod <- glm(group ~ x1 + x2 + x3 + x4 + x5 + x6, data=sim_dat,
                 family="binomial")

## Just check if function throws any errors
test_that("2 treatments, no conf_int, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aiptw",
                                           conf_int=F,
                                           outcome_model=outc_mod,
                                           treatment_model=treat_mod,
                                           cause=1), NA)
})

test_that("2 treatments, with conf_int, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aiptw",
                                           conf_int=T,
                                           outcome_model=outc_mod,
                                           treatment_model=treat_mod,
                                           cause=1), NA)
})

test_that("2 treatments, no conf_int, with boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aiptw",
                                           conf_int=F,
                                           bootstrap=T,
                                           n_boot=2,
                                           outcome_model=outc_mod,
                                           treatment_model=treat_mod,
                                           cause=1), NA)
})

test_that("2 treatments, no conf_int, no boot, with times", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aiptw",
                                           conf_int=F,
                                           bootstrap=F,
                                           n_boot=2,
                                           outcome_model=outc_mod,
                                           treatment_model=treat_mod,
                                           times=c(0.5, 1, 1.3),
                                           cause=1), NA)
})

cens_mod <- coxph(Surv(time, event==0) ~ x2, data=sim_dat, x=T)

test_that("2 treatments, no conf_int, no boot, with times, with cens_mod", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aiptw",
                                           conf_int=F,
                                           bootstrap=F,
                                           n_boot=2,
                                           outcome_model=outc_mod,
                                           treatment_model=treat_mod,
                                           censoring_model=cens_mod,
                                           times=c(0.5, 1, 1.3),
                                           cause=1), NA)
})