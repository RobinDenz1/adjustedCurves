library(adjustedCurves)
library(survival)
library(tmle)

set.seed(36)

sim_dat <- adjustedCurves::sim_confounded_surv(n=200, max_t=1.5)
sim_dat$time <- round(sim_dat$time * 10) + 1
sim_dat$event[sim_dat$event==1] <- sample(c(1, 2), size=sum(sim_dat$event),
                                          replace=T)

# outcome model
outcome_vars <- c("x1", "x2", "x3", "x4", "x5", "x6")

# treatment model
treat_mod <- glm(group ~ x2 + x3 + x5, data=sim_dat, family="binomial")

## Just check if function throws any errors
test_that("2 treatments, no conf_int, no boot, SL outcome, glm trt", {
  expect_error(suppressWarnings(adjustedcif(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="tmle_pseudo",
                                            conf_int=F,
                                            SL.ftime=c("SL.glm"),
                                            outcome_vars=outcome_vars,
                                            treatment_model=treat_mod,
                                            cause=1)), NA)
})

test_that("2 treatments, no conf_int, no boot, SL outcome, SL trt", {
  expect_error(suppressWarnings(adjustedcif(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="tmle_pseudo",
                                            conf_int=F,
                                            SL.ftime=c("SL.glm"),
                                            SL.trt=c("SL.glm"),
                                            treatment_vars=outcome_vars,
                                            outcome_vars=outcome_vars,
                                            treatment_model=NULL,
                                            cause=1)), NA)
})


test_that("2 treatments, with conf_int, no boot, SL outcome, glm trt", {
  expect_error(suppressWarnings(adjustedcif(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="tmle_pseudo",
                                            conf_int=T,
                                            SL.ftime=c("SL.glm"),
                                            outcome_vars=outcome_vars,
                                            treatment_model=treat_mod,
                                            cause=1)), NA)
})

test_that("2 treatments, with conf_int, no boot, SL outcome, SL trt", {
  expect_error(suppressWarnings(adjustedcif(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="tmle_pseudo",
                                            conf_int=T,
                                            SL.ftime=c("SL.glm"),
                                            SL.trt=c("SL.glm"),
                                            treatment_vars=outcome_vars,
                                            outcome_vars=outcome_vars,
                                            treatment_model=NULL,
                                            cause=1)), NA)
})


test_that("2 treatments, no conf_int, with boot, SL outcome, glm trt", {
  expect_error(suppressWarnings(adjustedcif(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="tmle_pseudo",
                                            conf_int=F,
                                            bootstrap=T,
                                            n_boot=2,
                                            SL.ftime=c("SL.glm"),
                                            outcome_vars=outcome_vars,
                                            treatment_model=treat_mod,
                                            cause=1)), NA)
})

test_that("2 treatments, no conf_int, with boot, SL outcome, SL trt", {
  expect_error(suppressWarnings(adjustedcif(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="tmle_pseudo",
                                            conf_int=F,
                                            bootstrap=T,
                                            n_boot=2,
                                            SL.ftime=c("SL.glm"),
                                            SL.trt=c("SL.glm"),
                                            treatment_vars=outcome_vars,
                                            outcome_vars=outcome_vars,
                                            treatment_model=NULL,
                                            cause=1)), NA)
})
