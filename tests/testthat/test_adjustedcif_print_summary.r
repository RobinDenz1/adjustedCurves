library(survival)

sim_dat <- sim_confounded_crisk(n=20)
sim_dat$group <- as.factor(sim_dat$group)

adj <- adjustedcif(data=sim_dat,
                   variable="group",
                   ev_time="time",
                   event="event",
                   method="aalen_johansen",
                   cause=1)

test_that("print.adjustedcif", {
  expect_error(print(adj), NA)
})

test_that("summary.adjustedcif, aalen_johansen", {
  expect_error(summary(adj), NA)
})

adj$method <- "direct"
test_that("summary.adjustedcif, direct", {
  expect_error(summary(adj), NA)
})

adj$method <- "direct_pseudo"
test_that("summary.adjustedcif, direct_pseudo", {
  expect_error(summary(adj), NA)
})

adj$method <- "iptw"
test_that("summary.adjustedcif, iptw", {
  expect_error(summary(adj), NA)
})

adj$method <- "iptw_pseudo"
test_that("summary.adjustedcif, iptw_pseudo", {
  expect_error(summary(adj), NA)
})

adj$method <- "matching"
test_that("summary.adjustedcif, matching", {
  expect_error(summary(adj), NA)
})

adj$method <- "aiptw"
test_that("summary.adjustedcif, aiptw", {
  expect_error(summary(adj), NA)
})

adj$method <- "aiptw_pseudo"
test_that("summary.adjustedcif, aiptw_pseudo", {
  expect_error(summary(adj), NA)
})

adj$method <- "tmle"
test_that("summary.adjustedcif, tmle", {
  expect_error(summary(adj), NA)
})

adj$call$bootstrap <- TRUE

test_that("summary.adjustedcif, with boot", {
  expect_error(summary(adj), NA)
})

adj$call$conf_int <- TRUE

test_that("summary.adjustedcif, with conf_int", {
  expect_error(summary(adj), NA)
})

adj$mids_analyses <- 1

test_that("summary.adjustedcif, with mids", {
  expect_error(summary(adj), NA)
})
