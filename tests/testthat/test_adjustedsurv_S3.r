
suppressMessages(requireNamespace("survival"))

set.seed(42)
sim_dat <- readRDS(system.file("testdata",
                               "d_sim_surv_n_20.Rds",
                               package="adjustedCurves"))
sim_dat$group <- as.factor(sim_dat$group)

adj <- adjustedsurv(data=sim_dat,
                    variable="group",
                    ev_time="time",
                    event="event",
                    method="km")

test_that("print.adjustedsurv", {
  expect_snapshot_output(print(adj))
})

test_that("summary.adjustedsurv, km", {
  expect_snapshot_output(summary(adj))
})

adj$method <- "direct"
test_that("summary.adjustedsurv, direct", {
  expect_snapshot_output(summary(adj))
})

adj$method <- "direct_pseudo"
test_that("summary.adjustedsurv, direct_pseudo", {
  expect_snapshot_output(summary(adj))
})

adj$method <- "iptw_km"
test_that("summary.adjustedsurv, iptw_km", {
  expect_snapshot_output(summary(adj))
})

adj$method <- "iptw_cox"
test_that("summary.adjustedsurv, iptw_cox", {
  expect_snapshot_output(summary(adj))
})

adj$method <- "iptw_pseudo"
test_that("summary.adjustedsurv, iptw_pseudo", {
  expect_snapshot_output(summary(adj))
})

adj$method <- "matching"
test_that("summary.adjustedsurv, matching", {
  expect_snapshot_output(summary(adj))
})

adj$method <- "emp_lik"
test_that("summary.adjustedsurv, emp_lik", {
  expect_snapshot_output(summary(adj))
})

adj$method <- "aiptw"
test_that("summary.adjustedsurv, aiptw", {
  expect_snapshot_output(summary(adj))
})

adj$method <- "aiptw_pseudo"
test_that("summary.adjustedsurv, aiptw_pseudo", {
  expect_snapshot_output(summary(adj))
})

adj$method <- "tmle"
test_that("summary.adjustedsurv, tmle", {
  expect_snapshot_output(summary(adj))
})

adj$method <- "strat_amato"
test_that("summary.adjustedsurv, strat_amato", {
  expect_snapshot_output(summary(adj))
})

adj$method <- "strat_nieto"
test_that("summary.adjustedsurv, strat_nieto", {
  expect_snapshot_output(summary(adj))
})

adj$method <- "strat_cupples"
test_that("summary.adjustedsurv, strat_cupples", {
  expect_snapshot_output(summary(adj))
})

adj$method <- "iv_2SRIF"
test_that("summary.adjustedsurv, iv_2SRIF", {
  expect_snapshot_output(summary(adj))
})

adj$method <- "prox_iptw"
test_that("summary.adjustedsurv, prox_iptw", {
  expect_snapshot_output(summary(adj))
})

adj$method <- "prox_aiptw"
test_that("summary.adjustedsurv, prox_aiptw", {
  expect_snapshot_output(summary(adj))
})

adj$call$bootstrap <- TRUE

test_that("summary.adjustedsurv, with boot", {
  expect_snapshot_output(summary(adj))
})

adj$call$conf_int <- TRUE

test_that("summary.adjustedsurv, with conf_int", {
  expect_snapshot_output(summary(adj))
})

adj$mids_analyses <- 1

test_that("summary.adjustedsurv, with mids", {
  expect_snapshot_output(summary(adj))
})
