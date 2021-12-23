
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
  adj_test <- adjusted_curve_diff(adj, to=1.3)
  expect_equal(round(adj_test$observed_diff_integral, 4), 0.2197)
  expect_equal(round(adj_test$integral_se, 4), 0.0474)
  expect_equal(round(adj_test$p_value, 4), 0)
  expect_equal(adj_test$n_boot, 10)
})

test_that("survival, 2 treatments, with from", {
  adj_test <- adjusted_curve_diff(adj, to=1.3, from=0.5)
  expect_equal(round(adj_test$observed_diff_integral, 4), 0.1682)
  expect_equal(round(adj_test$integral_se, 4), 0.0394)
  expect_equal(round(adj_test$p_value, 4), 0)
  expect_equal(adj_test$n_boot, 10)
})

sim_dat$group <- as.character(sim_dat$group)
sim_dat$group[sim_dat$group=="1"] <- sample(x=c(1, 2),
                                  size=nrow(sim_dat[sim_dat$group=="1", ]),
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
  adj_test <- adjusted_curve_diff(adj, to=1.3)
  expect_equal(round(adj_test$`0 vs. 1`$observed_diff_integral, 4), 0.1465)
  expect_equal(round(adj_test$`0 vs. 1`$integral_se, 4), 0.1021)
  expect_equal(round(adj_test$`0 vs. 1`$p_value, 4), 0.1111)
  expect_equal(adj_test$`0 vs. 1`$n_boot, 9)
})

test_that("survival, > 2 treatments, with from", {
  adj_test <- adjusted_curve_diff(adj, to=1.3, from=0.5)
  expect_equal(round(adj_test$`0 vs. 1`$observed_diff_integral, 4), 0.1126)
  expect_equal(round(adj_test$`0 vs. 1`$integral_se, 4), 0.0809)
  expect_equal(round(adj_test$`0 vs. 1`$p_value, 4), 0.1111)
  expect_equal(adj_test$`0 vs. 1`$n_boot, 9)
})

# competing risks case

sim_dat <- sim_confounded_crisk(n=200, max_t=1.5)
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
  adj_test <- adjusted_curve_diff(adj, to=1)
  expect_equal(round(adj_test$observed_diff_integral, 4), 0.1378)
  expect_equal(round(adj_test$integral_se, 4), 0.116)
  expect_equal(round(adj_test$p_value, 4), 0.2)
  expect_equal(adj_test$n_boot, 10)
})

test_that("CIF, 2 treatments, with from", {
  adj_test <- adjusted_curve_diff(adj, to=1, from=0.5)
  expect_equal(round(adj_test$observed_diff_integral, 4), 0.0816)
  expect_equal(round(adj_test$integral_se, 4), 0.068)
  expect_equal(round(adj_test$p_value, 4), 0.2)
  expect_equal(adj_test$n_boot, 10)
})

sim_dat$group <- as.character(sim_dat$group)
sim_dat$group[sim_dat$group=="1"] <- sample(x=c(1, 2),
                                    size=nrow(sim_dat[sim_dat$group=="1", ]),
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
  adj_test <- adjusted_curve_diff(adj, to=1)
  expect_equal(round(adj_test$`1 vs. 0`$observed_diff_integral, 4), 0.1386)
  expect_equal(round(adj_test$`1 vs. 0`$integral_se, 4), 0.1427)
  expect_equal(round(adj_test$`1 vs. 0`$p_value, 4), 0.3)
  expect_equal(adj_test$`1 vs. 0`$n_boot, 10)
})

test_that("CIF, > 2 treatments, with from", {
  adj_test <- adjusted_curve_diff(adj, to=1, from=0.5)
  expect_equal(round(adj_test$`1 vs. 0`$observed_diff_integral, 4), 0.0793)
  expect_equal(round(adj_test$`1 vs. 0`$integral_se, 4), 0.0838)
  expect_equal(round(adj_test$`1 vs. 0`$p_value, 4), 0.3)
  expect_equal(adj_test$`1 vs. 0`$n_boot, 10)
})
