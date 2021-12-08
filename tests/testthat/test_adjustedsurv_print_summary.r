library(survival)

sim_dat <- sim_confounded_surv(n=20)
sim_dat$group <- as.factor(sim_dat$group)

adj <- adjustedsurv(data=sim_dat,
                    variable="group",
                    ev_time="time",
                    event="event",
                    method="km")

test_that("print.adjustedsurv", {
  expect_error(print(adj), NA)
})

test_that("summary.adjustedsurv, km", {
  expect_error(summary(adj), NA)
})

adj$method <- "direct"
test_that("summary.adjustedsurv, direct", {
  expect_error(summary(adj), NA)
})

adj$method <- "direct_pseudo"
test_that("summary.adjustedsurv, direct_pseudo", {
  expect_error(summary(adj), NA)
})

adj$method <- "iptw_km"
test_that("summary.adjustedsurv, iptw_km", {
  expect_error(summary(adj), NA)
})

adj$method <- "iptw_cox"
test_that("summary.adjustedsurv, iptw_cox", {
  expect_error(summary(adj), NA)
})

adj$method <- "iptw_pseudo"
test_that("summary.adjustedsurv, iptw_pseudo", {
  expect_error(summary(adj), NA)
})

adj$method <- "matching"
test_that("summary.adjustedsurv, matching", {
  expect_error(summary(adj), NA)
})

adj$method <- "emp_lik"
test_that("summary.adjustedsurv, emp_lik", {
  expect_error(summary(adj), NA)
})

adj$method <- "aiptw"
test_that("summary.adjustedsurv, aiptw", {
  expect_error(summary(adj), NA)
})

adj$method <- "aiptw_pseudo"
test_that("summary.adjustedsurv, aiptw_pseudo", {
  expect_error(summary(adj), NA)
})

adj$method <- "tmle"
test_that("summary.adjustedsurv, tmle", {
  expect_error(summary(adj), NA)
})

adj$method <- "ostmle"
test_that("summary.adjustedsurv, ostmle", {
  expect_error(summary(adj), NA)
})

adj$method <- "strat_amato"
test_that("summary.adjustedsurv, strat_amato", {
  expect_error(summary(adj), NA)
})

adj$method <- "strat_gregory_nieto"
test_that("summary.adjustedsurv, strat_gregory_nieto", {
  expect_error(summary(adj), NA)
})

adj$method <- "strat_cupples"
test_that("summary.adjustedsurv, strat_cupples", {
  expect_error(summary(adj), NA)
})

adj$call$bootstrap <- TRUE

test_that("summary.adjustedsurv, with boot", {
  expect_error(summary(adj), NA)
})

adj$call$conf_int <- TRUE

test_that("summary.adjustedsurv, with conf_int", {
  expect_error(summary(adj), NA)
})

adj$mids_analyses <- 1

test_that("summary.adjustedsurv, with mids", {
  expect_error(summary(adj), NA)
})
