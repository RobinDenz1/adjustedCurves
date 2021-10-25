library(survival)
library(cmprsk)
library(adjustedCurves)

# survival case

sim_dat <- sim_confounded_surv(n=200, max_t=1.5)
sim_dat$group <- factor(sim_dat$group)

adj <- adjustedsurv(data=sim_dat,
                    variable="group",
                    ev_time="time",
                    event="event",
                    method="km",
                    conf_int=TRUE,
                    bootstrap=TRUE,
                    n_boot=10)

test_that("test_curve_equality, no from", {
  expect_error(test_curve_equality(adj, to=1.3), NA)
})

test_that("test_curve_equality, with from", {
  expect_error(test_curve_equality(adj, to=1.3, from=0.5), NA)
})

# competing risks case

sim_dat <- sim_confounded_surv(n=200, max_t=1.5)
sim_dat$event[sim_dat$event==1] <- sample(c(1, 2), size=sum(sim_dat$event),
                                          replace=TRUE)
sim_dat$group <- factor(sim_dat$group)

adj <- adjustedcif(data=sim_dat,
                   variable="group",
                   ev_time="time",
                   event="event",
                   method="aalen_johansen",
                   cause=1,
                   conf_int=TRUE,
                   bootstrap=TRUE,
                   n_boot=10)

test_that("test_curve_equality, no from", {
  expect_error(test_curve_equality(adj, to=1.3), NA)
})

test_that("test_curve_equality, with from", {
  expect_error(test_curve_equality(adj, to=1.3, from=0.5), NA)
})
