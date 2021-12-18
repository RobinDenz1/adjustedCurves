library(dplyr)

set.seed(42)

sim_dat <- adjustedCurves::sim_confounded_surv(n=50)
sim_dat$group <- as.factor(sim_dat$group)

sim_dat_tibble <- tibble(sim_dat)

# adjustedsurv works with tibbles?
test_that("adjustedsurv, tibbles, no boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat_tibble,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="iptw_km",
                                            conf_int=FALSE,
                                            treatment_model=group ~ x1), NA)
})

test_that("adjustedsurv, tibbles, with boot", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat_tibble,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="iptw_km",
                                            bootstrap=TRUE,
                                            n_boot=2,
                                            conf_int=FALSE,
                                            treatment_model=group ~ x1), NA)
})

# adjustedcif works with tibbles?
test_that("adjustedcif, tibbles, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat_tibble,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aalen_johansen",
                                           conf_int=FALSE,
                                           cause=1), NA)
})

test_that("adjustedcif, tibbles, with boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat_tibble,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aalen_johansen",
                                           bootstrap=TRUE,
                                           n_boot=2,
                                           conf_int=FALSE,
                                           cause=1), NA)
})
