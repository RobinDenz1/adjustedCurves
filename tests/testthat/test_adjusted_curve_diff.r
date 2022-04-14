
suppressMessages(requireNamespace("survival"))

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
                    bootstrap=FALSE)

sim_dat$event[1] <- 2
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="aalen_johansen",
                      conf_int=TRUE,
                      bootstrap=TRUE,
                      n_boot=2,
                      cause=1)

test_that("surv, default", {
  adj_diff <- adjusted_curve_diff(adj)
  expect_equal(nrow(adj_diff), 193)
  expect_equal(round(adj_diff$diff[1], 4), 0)
  expect_equal(round(adj_diff$diff[139], 4), -0.291)
})

test_that("surv, times", {
  adj_diff <- adjusted_curve_diff(adj, times=0.8)
  expect_equal(nrow(adj_diff), 1)
  expect_equal(round(adj_diff$diff[1], 4), -0.2425)
})

test_that("cif, times", {
  adj_diff <- adjusted_curve_diff(adjcif, times=0.8, conf_int=TRUE,
                                  use_boot=TRUE)
  expect_equal(nrow(adj_diff), 1)
  expect_equal(round(adj_diff$diff[1], 4), 0.2343)
})
