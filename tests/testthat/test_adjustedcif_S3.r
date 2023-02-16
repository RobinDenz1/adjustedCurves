
suppressMessages(requireNamespace("survival"))

set.seed(42)
sim_dat <- readRDS(system.file("testdata",
                               "d_sim_crisk_n_20.Rds",
                               package="adjustedCurves"))
sim_dat$group <- as.factor(sim_dat$group)

adj <- adjustedcif(data=sim_dat,
                   variable="group",
                   ev_time="time",
                   event="event",
                   method="aalen_johansen",
                   cause=1)

test_that("print.adjustedcif", {
  expect_snapshot_output(print(adj))
})

test_that("summary.adjustedcif, aalen_johansen", {
  expect_snapshot_output(summary(adj))
})

adj$method <- "direct"
test_that("summary.adjustedcif, direct", {
  expect_snapshot_output(summary(adj))
})

adj$method <- "direct_pseudo"
test_that("summary.adjustedcif, direct_pseudo", {
  expect_snapshot_output(summary(adj))
})

adj$method <- "iptw"
test_that("summary.adjustedcif, iptw", {
  expect_snapshot_output(summary(adj))
})

adj$method <- "iptw_pseudo"
test_that("summary.adjustedcif, iptw_pseudo", {
  expect_snapshot_output(summary(adj))
})

adj$method <- "matching"
test_that("summary.adjustedcif, matching", {
  expect_snapshot_output(summary(adj))
})

adj$method <- "aiptw"
test_that("summary.adjustedcif, aiptw", {
  expect_snapshot_output(summary(adj))
})

adj$method <- "aiptw_pseudo"
test_that("summary.adjustedcif, aiptw_pseudo", {
  expect_snapshot_output(summary(adj))
})

adj$call$bootstrap <- TRUE

test_that("summary.adjustedcif, with boot", {
  expect_snapshot_output(summary(adj))
})

adj$call$conf_int <- TRUE

test_that("summary.adjustedcif, with conf_int", {
  expect_snapshot_output(summary(adj))
})

adj$mids_analyses <- 1

test_that("summary.adjustedcif, with mids", {
  expect_snapshot_output(summary(adj))
})
