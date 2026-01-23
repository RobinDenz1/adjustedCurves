
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
                    bootstrap=TRUE,
                    n_boot=10)

test_that("survival, 2 treatments, no from", {
  adj_test <- adjusted_curve_test(adj, to=1.3)
  expect_equal(round(adj_test$observed_diff_integral, 4), -0.2197)
  expect_equal(round(adj_test$integral_se, 4), 0.0474)
  expect_equal(round(adj_test$p_value, 4), 0)
  expect_equal(adj_test$n_boot, 10)
})

test_that("survival, 2 treatments, no from, linear", {
  adj_test <- adjusted_curve_test(adj, to=1.3, interpolation="linear")
  expect_equal(round(adj_test$observed_diff_integral, 4), -0.2193)
  expect_equal(round(adj_test$integral_se, 4), 0.0469)
  expect_equal(round(adj_test$p_value, 4), 0)
  expect_equal(adj_test$n_boot, 10)
})

test_that("survival, 2 treatments, with from", {
  adj_test <- adjusted_curve_test(adj, to=1.3, from=0.5)
  expect_equal(round(adj_test$observed_diff_integral, 4), -0.1682)
  expect_equal(round(adj_test$integral_se, 4), 0.0394)
  expect_equal(round(adj_test$p_value, 4), 0)
  expect_equal(adj_test$n_boot, 10)
})

test_that("survival, 2 treatments, S3 methods", {
  adj_test <- adjusted_curve_test(adj, to=1.3, from=0.5)
  expect_snapshot_output(print(adj_test))
  expect_snapshot_output(summary(adj_test))
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
  adj_test <- adjusted_curve_test(adj, to=1.3)
  expect_equal(round(adj_test$`0 vs. 1`$observed_diff_integral, 4), -0.1465)
  expect_equal(round(adj_test$`0 vs. 1`$integral_se, 4), 0.0982)
  expect_equal(round(adj_test$`0 vs. 1`$p_value, 4), 0.1)
  expect_equal(adj_test$`0 vs. 1`$n_boot, 10)
})

test_that("survival, > 2 treatments, no from, linear", {
  adj_test <- adjusted_curve_test(adj, to=1.3, interpolation="linear")
  expect_equal(round(adj_test$`0 vs. 1`$observed_diff_integral, 4), -0.1467)
  expect_equal(round(adj_test$`0 vs. 1`$integral_se, 4), 0.0981)
  expect_equal(round(adj_test$`0 vs. 1`$p_value, 4), 0.1)
  expect_equal(adj_test$`0 vs. 1`$n_boot, 10)
})

test_that("survival, > 2 treatments, with from", {
  adj_test <- adjusted_curve_test(adj, to=1.3, from=0.5)
  expect_equal(round(adj_test$`0 vs. 1`$observed_diff_integral, 4), -0.1126)
  expect_equal(round(adj_test$`0 vs. 1`$integral_se, 4), 0.0764)
  expect_equal(round(adj_test$`0 vs. 1`$p_value, 4), 0.10)
  expect_equal(adj_test$`0 vs. 1`$n_boot, 10)
})

test_that("survival, > 2 treatments, S3 method", {
  adj_test <- adjusted_curve_test(adj, to=1.3, from=0.5)
  adj_test$method <- "iptw_km"
  expect_snapshot_output(print(adj_test))
  expect_snapshot_output(summary(adj_test))
})

# competing risks case

set.seed(1234)
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
  set.seed(134)
  adj_test <- adjusted_curve_test(adj, to=1)
  expect_equal(round(adj_test$observed_diff_integral, 4), -0.0164)
  expect_equal(round(adj_test$integral_se, 4), 0.0453)
  expect_equal(round(adj_test$p_value, 4), 0.8)
  expect_equal(adj_test$n_boot, 10)
})

test_that("CIF, 2 treatments, no from, linear", {
  set.seed(134)
  adj_test <- adjusted_curve_test(adj, to=1, interpolation="linear")
  expect_equal(round(adj_test$observed_diff_integral, 4), -0.0154)
  expect_equal(round(adj_test$integral_se, 4), 0.0458)
  expect_equal(round(adj_test$p_value, 4), 0.8)
  expect_equal(adj_test$n_boot, 10)
})

test_that("CIF, 2 treatments, with from", {
  set.seed(1234)
  adj_test <- adjusted_curve_test(adj, to=1, from=0.5)
  expect_equal(round(adj_test$observed_diff_integral, 3), -0.001)
  expect_equal(round(adj_test$integral_se, 4), 0.0296)
  expect_equal(round(adj_test$p_value, 4), 1)
  expect_equal(adj_test$n_boot, 10)
})

test_that("CIF, 2 treatments, S3 methods", {
  set.seed(455)
  adj_test <- adjusted_curve_test(adj, to=1.3, from=0.5)
  adj_test$method <- "iptw"
  expect_snapshot_output(print(adj_test))
  expect_snapshot_output(summary(adj_test))
})

set.seed(1324)
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
  set.seed(1234)
  adj_test <- adjusted_curve_test(adj, to=1)
  expect_equal(round(adj_test$`0 vs. 1`$observed_diff_integral, 3), 0)
  expect_equal(round(adj_test$`0 vs. 2`$integral_se, 4), 0.0539)
  expect_equal(round(adj_test$`1 vs. 2`$p_value, 4), 0.7)
  expect_equal(adj_test$`0 vs. 1`$n_boot, 10)
})

test_that("CIF, > 2 treatments, no from, linear", {
  set.seed(4525)
  adj_test <- adjusted_curve_test(adj, to=1, interpolation="linear")
  expect_equal(round(adj_test$`0 vs. 1`$observed_diff_integral, 4), 0.0016)
  expect_equal(round(adj_test$`0 vs. 2`$integral_se, 4), 0.054)
  expect_equal(round(adj_test$`1 vs. 2`$p_value, 4), 0.7)
  expect_equal(adj_test$`0 vs. 1`$n_boot, 10)
})

test_that("CIF, > 2 treatments, with from", {
  set.seed(1234)
  adj_test <- adjusted_curve_test(adj, to=1, from=0.5)
  expect_equal(round(adj_test$`0 vs. 1`$observed_diff_integral, 4), 0.0073)
  expect_equal(round(adj_test$`0 vs. 2`$integral_se, 4), 0.0308)
  expect_equal(round(adj_test$`1 vs. 2`$p_value, 4), 0.7)
  expect_equal(adj_test$`0 vs. 1`$n_boot, 10)
})

test_that("CIF, > 2 treatments, S3 methods", {
  set.seed(12345)
  adj_test <- adjusted_curve_test(adj, to=1.3, from=0.5)
  expect_snapshot_output(print(adj_test))
  expect_snapshot_output(summary(adj_test))
})
