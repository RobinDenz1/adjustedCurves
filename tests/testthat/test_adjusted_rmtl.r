library(survival)
library(riskRegression)
library(prodlim)

### using single-event data

set.seed(42)
sim_dat <- sim_confounded_surv(n=50, max_t=1.5)
sim_dat$group <- factor(sim_dat$group)

adj <- adjustedsurv(data=sim_dat,
                    variable="group",
                    ev_time="time",
                    event="event",
                    method="km",
                    conf_int=TRUE,
                    bootstrap=TRUE,
                    n_boot=10)

test_that("rmtl surv, no boot", {
  expect_error(adjusted_rmtl(adj, to=1.1), NA)
})

test_that("rmtl surv, with boot", {
  expect_error(adjusted_rmtl(adj, to=1.1, use_boot=TRUE), NA)
})

test_that("rmtl surv, no boot, using from", {
  expect_error(adjusted_rmtl(adj, to=1.1, from=0.3), NA)
})

test_that("rmtl surv, with boot, using from", {
  expect_error(adjusted_rmtl(adj, to=1.1, from=0.3, use_boot=TRUE), NA)
})

### using competing risks data

set.seed(42)
sim_dat <- sim_confounded_crisk(n=50, max_t=1.5)
sim_dat$group <- factor(sim_dat$group)

adj <- adjustedcif(data=sim_dat,
                   variable="group",
                   ev_time="time",
                   event="event",
                   method="aalen_johansen",
                   conf_int=TRUE,
                   bootstrap=TRUE,
                   n_boot=10,
                   cause=1)

test_that("rmtl cif, no boot", {
  expect_error(adjusted_rmtl(adj, to=1.1), NA)
})

test_that("rmtl cif, with boot", {
  expect_error(adjusted_rmtl(adj, to=1.1, use_boot=TRUE), NA)
})

test_that("rmtl cif, no boot, using from", {
  expect_error(adjusted_rmtl(adj, to=1.1, from=0.3), NA)
})

test_that("rmtl cif, with boot, using from", {
  expect_error(adjusted_rmtl(adj, to=1.1, from=0.3, use_boot=TRUE), NA)
})
