
suppressMessages(requireNamespace("survival"))

set.seed(42)
sim_dat <- readRDS(system.file("testdata",
                               "d_sim_surv_n_50.Rds",
                               package="adjustedCurves"))
sim_dat$group <- factor(sim_dat$group)

adj <- adjustedsurv(data=sim_dat,
                    variable="group",
                    ev_time="time",
                    event="event",
                    method="km",
                    conf_int=TRUE,
                    bootstrap=TRUE,
                    n_boot=2)

test_that("q surv defaults", {
  adj_q <- adjusted_surv_quantile(adj)
  expect_equal(round(adj_q$q_surv, 4), c(0.4785, 0.6252))
})

test_that("q surv, conf_int, no boot", {
  adj_q <- adjusted_surv_quantile(adj, conf_int=TRUE)
  expect_equal(round(adj_q$q_surv, 4), c(0.4785, 0.6252))
  expect_equal(round(adj_q$ci_lower, 4), c(0.3848, 0.5115))
  expect_equal(round(adj_q$ci_upper, 4), c(0.9498, 1.0331))
})

test_that("q surv, conf_int, with boot", {
  adj_q <- adjusted_surv_quantile(adj, conf_int=TRUE, use_boot=TRUE)
  expect_equal(round(adj_q$q_surv, 4), c(0.4785, 0.6252))
  expect_equal(round(adj_q$ci_lower, 4), c(0.4020, 0.7408))
  expect_equal(round(adj_q$ci_upper, 4), c(0.4785, 0.7502))
})

test_that("q surv, conf_int, with boot, linear", {
  adj_q <- adjusted_surv_quantile(adj, conf_int=TRUE, use_boot=TRUE,
                                  interpolation="linear")
  expect_equal(round(adj_q$q_surv, 4), c(0.4832, 0.6279))
  expect_equal(round(adj_q$ci_lower, 4), c(0.4204, 0.7165))
  expect_equal(round(adj_q$ci_upper, 4), c(0.4971, 0.7437))
})

test_that("q surv, multiple p", {
  adj_q <- adjusted_surv_quantile(adj, p=c(0.1, 0.2), conf_int=TRUE)
  expect_equal(round(adj_q$q_surv, 4), c(1.0503, 0.6561, 1.2304, 0.8275))
  expect_equal(round(adj_q$ci_lower, 4), c(0.6561, 0.6208, 0.8275, 0.7408))
  expect_equal(sum(is.na(adj_q$ci_upper)), 4)
})

