library(survival)
library(adjustedCurves)

set.seed(42)

sim_dat <- sim_confounded_surv(n=20, max_t=1.5)
sim_dat$group <- factor(sim_dat$group)

adj <- adjustedsurv(data=sim_dat,
                    variable="group",
                    ev_time="time",
                    event="event",
                    method="km",
                    conf_int=TRUE,
                    bootstrap=TRUE,
                    n_boot=2)

adj_test <- adjusted_curve_diff(adj, to=0.2)

test_that("print.curve_test, default", {
  expect_error(print(adj_test), NA)
})

test_that("summary.curve_test, default", {
  expect_error(summary(adj_test), NA)
})
