library(nnet)

set.seed(42)

sim_dat <- adjustedCurves::sim_confounded_surv(n=50)
sim_dat$group <- as.factor(sim_dat$group)

mod <- glm(group ~ x1 + x2 + x3 + x4 + x5 + x6, data=sim_dat,
           family="binomial")

## Just check if function throws any errors
test_that("2 treatments, no conf_int, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="iptw_km",
                                            conf_int=F,
                                            treatment_model=mod), NA)
})

test_that("2 treatments, with conf_int, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="iptw_km",
                                            conf_int=T,
                                            treatment_model=mod), NA)
})

test_that("2 treatments, no conf_int, with boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="iptw_km",
                                            conf_int=F,
                                            bootstrap=T,
                                            n_boot=2,
                                            treatment_model=mod), NA)
})

test_that("2 treatments, with conf_int, with boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="iptw_km",
                                            conf_int=T,
                                            bootstrap=T,
                                            n_boot=2,
                                            treatment_model=mod), NA)
})

test_that("2 treatments, no conf_int, with WeightIt", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="iptw_km",
                                            conf_int=F,
                                            bootstrap=F,
                                            treatment_model=group ~ x1 + x2,
                                            weight_method="ps"), NA)
})

test_that("2 treatments, no conf_int, with user-weights", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="iptw_km",
                                            conf_int=F,
                                            bootstrap=F,
                                            treatment_model=runif(n=50, min=1, max=2))
               , NA)
})

sim_dat <- adjustedCurves::sim_confounded_surv(n=50)
sim_dat$group[sim_dat$group==1] <- sample(c(1, 2),
                                          size=nrow(sim_dat[sim_dat$group==1,]),
                                          replace=T)
sim_dat$group <- as.factor(sim_dat$group)

mod <- nnet::multinom(group ~ x1 + x2 + x3 + x4 + x5 + x6, data=sim_dat)

test_that("> 2 treatments, no conf_int, no boot, no ...", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="iptw_km",
                                            conf_int=F,
                                            treatment_model=mod), NA)
})

test_that("> 2 treatments, with conf_int, no boot, no ...", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="iptw_km",
                                            conf_int=T,
                                            treatment_model=mod), NA)
})

test_that("> 2 treatments, no conf_int, with boot, no ...", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="iptw_km",
                                            conf_int=F,
                                            bootstrap=T,
                                            n_boot=2,
                                            treatment_model=mod), NA)
})

test_that("> 2 treatments, with conf_int, with boot, no ...", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="iptw_km",
                                            conf_int=T,
                                            bootstrap=T,
                                            n_boot=2,
                                            treatment_model=mod), NA)
})

# DONT RUN: would require package dependency on "mlogit"
#test_that("> 2 treatments, no conf_int, with WeightIt", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="iptw_km",
#                                            conf_int=F,
#                                            bootstrap=F,
#                                            treatment_model=group ~ x1 + x2,
#                                            weight_method="ps"), NA)
#})

test_that("> 2 treatments, no conf_int, with user-weights", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="iptw_km",
                                            conf_int=F,
                                            bootstrap=F,
                                            treatment_model=runif(n=50, min=1, max=2))
               , NA)
})

