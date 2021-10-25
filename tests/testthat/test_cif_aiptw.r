library(survival)
library(riskRegression)

set.seed(42)

sim_dat <- adjustedCurves::sim_confounded_surv(n=50)
sim_dat$event[sim_dat$event==1] <- sample(c(1, 2), size=sum(sim_dat$event),
                                          replace=TRUE)
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
                                           conf_int=FALSE,
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
                                           conf_int=TRUE,
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
                                           conf_int=FALSE,
                                           bootstrap=TRUE,
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
                                           conf_int=FALSE,
                                           bootstrap=FALSE,
                                           n_boot=2,
                                           outcome_model=outc_mod,
                                           treatment_model=treat_mod,
                                           times=c(0.5, 1),
                                           cause=1), NA)
})

cens_mod <- coxph(Surv(time, event==0) ~ x2, data=sim_dat, x=TRUE)

test_that("2 treatments, no conf_int, no boot, with times, with cens_mod", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aiptw",
                                           conf_int=FALSE,
                                           bootstrap=FALSE,
                                           n_boot=2,
                                           outcome_model=outc_mod,
                                           treatment_model=treat_mod,
                                           censoring_model=cens_mod,
                                           times=c(0.5, 1),
                                           cause=1), NA)
})
