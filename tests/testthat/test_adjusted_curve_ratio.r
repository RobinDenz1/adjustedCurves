
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

# competing risks case
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
  adj_diff <- adjusted_curve_ratio(adj)
  expect_equal(nrow(adj_diff), 193)
  expect_equal(round(adj_diff$ratio[1], 4), 1)
  expect_equal(round(adj_diff$ratio[139], 4), 0.5049)
})

test_that("surv, times", {
  adj_diff <- adjusted_curve_ratio(adj, times=0.8)
  expect_equal(nrow(adj_diff), 1)
  expect_equal(round(adj_diff$ratio[1], 4), 0.4548)
})

test_that("surv, times, conf_int", {
  adj_diff <- adjusted_curve_ratio(adj, times=0.8, conf_int=TRUE)
  expect_equal(nrow(adj_diff), 1)
  expect_equal(round(adj_diff$ratio[1], 4), 0.4548)
  expect_equal(round(adj_diff$ci_lower[1], 3), -0.585)
  expect_equal(round(adj_diff$ci_upper[1], 3), 1.999)
  expect_equal(round(adj_diff$p_value[1], 3), 0.352)
})

test_that("cif, times", {
  adj_diff <- adjusted_curve_ratio(adjcif, times=0.8)
  expect_equal(nrow(adj_diff), 1)
  expect_equal(round(adj_diff$ratio[1], 4), 1.4344)
})

test_that("cif, times, conf_int", {
  adj_diff <- adjusted_curve_ratio(adjcif, times=0.8, conf_int=TRUE,
                                   use_boot=TRUE)
  expect_equal(nrow(adj_diff), 1)
  expect_equal(round(adj_diff$ratio[1], 4), 1.4344)
  expect_equal(round(adj_diff$ci_lower[1], 4), 1.3341)
  expect_equal(round(adj_diff$ci_upper[1], 4), 1.5466)
  expect_equal(round(adj_diff$p_value[1], 4), 0)
})
