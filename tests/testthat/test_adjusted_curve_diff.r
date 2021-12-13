library(survival)
library(cmprsk)
library(adjustedCurves)

# survival case

set.seed(42)

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

test_that("survival, 2 treatments, no from", {
  expect_error(adjusted_curve_diff(adj, to=1.3), NA)
})

test_that("survival, 2 treatments, with from", {
  expect_error(adjusted_curve_diff(adj, to=1.3, from=0.5), NA)
})

sim_dat$group <- as.character(sim_dat$group)
sim_dat$group[sim_dat$group=="1"] <- sample(x=c(1, 2),
                                  size=nrow(sim_dat[sim_dat$group=="1",]),
                                  replace=TRUE)
sim_dat$group <- as.factor(sim_dat$group)

adj <- adjustedsurv(data=sim_dat,
                    variable="group",
                    ev_time="time",
                    event="event",
                    method="km",
                    conf_int=TRUE,
                    bootstrap=TRUE,
                    n_boot=10)

test_that("survival, > 2 treatments, no from", {
  expect_error(adjusted_curve_diff(adj, to=1.3), NA)
})

test_that("survival, > 2 treatments, with from", {
  expect_error(adjusted_curve_diff(adj, to=1.3, from=0.5), NA)
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

test_that("CIF, 2 treatments, no from", {
  expect_error(adjusted_curve_diff(adj, to=1), NA)
})

test_that("CIF, 2 treatments, with from", {
  expect_error(adjusted_curve_diff(adj, to=1, from=0.5), NA)
})

sim_dat$group <- as.character(sim_dat$group)
sim_dat$group[sim_dat$group=="1"] <- sample(x=c(1, 2),
                                            size=nrow(sim_dat[sim_dat$group=="1",]),
                                            replace=TRUE)
sim_dat$group <- as.factor(sim_dat$group)

adj <- adjustedcif(data=sim_dat,
                   variable="group",
                   ev_time="time",
                   event="event",
                   method="aalen_johansen",
                   cause=1,
                   conf_int=TRUE,
                   bootstrap=TRUE,
                   n_boot=10)

test_that("CIF, > 2 treatments, no from", {
  expect_error(adjusted_curve_diff(adj, to=1), NA)
})

test_that("CIF, > 2 treatments, with from", {
  expect_error(adjusted_curve_diff(adj, to=1, from=0.5), NA)
})

