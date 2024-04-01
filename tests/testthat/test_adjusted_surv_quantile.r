
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

test_that("q surv defaults, difference", {
  adj_q <- adjusted_surv_quantile(adj, contrast="diff")
  expect_equal(round(adj_q$diff, 4), -0.1468)
})

test_that("q surv defaults, difference, groups", {
  adj_q <- adjusted_surv_quantile(adj, contrast="diff", group_1="1",
                                  group_2="0")
  expect_equal(round(adj_q$diff, 4), 0.1468)
})

test_that("q surv defaults, ratio", {
  adj_q <- adjusted_surv_quantile(adj, contrast="ratio")
  expect_equal(round(adj_q$ratio, 4), 0.7653)
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

test_that("q surv, conf_int, difference", {
  adj_q <- adjusted_surv_quantile(adj, conf_int=TRUE, contrast="diff")
  expect_equal(round(adj_q$diff, 4), -0.1468)
  expect_equal(round(adj_q$se, 4), 0.1165)
  expect_equal(round(adj_q$ci_lower, 4), -0.3751)
  expect_equal(round(adj_q$ci_upper, 4), 0.0816)
  expect_equal(round(adj_q$p_value, 4), 0.2077)
})

test_that("q surv, conf_int, difference, group", {
  adj_q <- adjusted_surv_quantile(adj, conf_int=TRUE, contrast="diff",
                                  group_1="1", group_2="0")
  expect_equal(round(adj_q$diff, 4), 0.1468)
  expect_equal(round(adj_q$ci_lower, 4), -0.0816)
  expect_equal(round(adj_q$ci_upper, 4), 0.3751)
})

test_that("q surv, conf_int, ratio", {
  adj_q <- adjusted_surv_quantile(adj, conf_int=TRUE, contrast="ratio")
  expect_equal(round(adj_q$ratio, 4), 0.7653)
  expect_equal(round(adj_q$ci_lower, 4), 0.4006)
  expect_equal(round(adj_q$ci_upper, 4), 1.1306)
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

test_that("q surv, multiple p, difference", {
  adj_q <- adjusted_surv_quantile(adj, p=c(0.1, 0.2), conf_int=TRUE,
                                  contrast="diff")
  expect_equal(round(adj_q$diff, 4), c(-0.1801, -0.1714))
  expect_equal(round(adj_q$ci_lower, 4), c(-0.1801, -0.6733))
  expect_equal(round(adj_q$ci_upper, 4), c(-0.1801, 0.3306))
})

test_that("q surv, multiple p, ratio", {
  adj_q <- adjusted_surv_quantile(adj, p=c(0.1, 0.2), conf_int=TRUE,
                                  contrast="ratio")
  expect_equal(round(adj_q$ratio, 4), c(0.8536, 0.7929))
  expect_equal(round(adj_q$ci_lower, 4), c(0.8536,  0.2761))
  expect_equal(round(adj_q$ci_upper, 4), c(0.8536, 1.5041))
})

## changing confidence levels is only supported for methods other than "km"

ps_mod <- glm(group ~ x1 + x4, data=sim_dat, family="binomial")

set.seed(4355)

adj <- adjustedsurv(data=sim_dat,
                    variable="group",
                    ev_time="time",
                    event="event",
                    method="iptw_km",
                    treatment_model=ps_mod,
                    conf_int=TRUE,
                    bootstrap=TRUE,
                    n_boot=2)

set.seed(4355)

adj0.9 <- adjustedsurv(data=sim_dat,
                       variable="group",
                       ev_time="time",
                       event="event",
                       method="iptw_km",
                       treatment_model=ps_mod,
                       conf_int=TRUE,
                       conf_level=0.9,
                       bootstrap=TRUE,
                       n_boot=2)

test_that("q surv, conf_int, no boot, conf_level", {
  adj_q <- adjusted_surv_quantile(adj, conf_int=TRUE, conf_level=0.9)
  adj_q0.9 <- adjusted_surv_quantile(adj0.9, conf_int=TRUE)
  expect_equal(adj_q$q_surv, adj_q0.9$q_surv)
  expect_equal(adj_q$ci_lower, adj_q0.9$ci_lower)
  expect_equal(adj_q$ci_upper, adj_q0.9$ci_upper)
})

test_that("q surv, conf_int, with boot, conf_level", {
  adj_q <- adjusted_surv_quantile(adj, conf_int=TRUE, conf_level=0.9,
                                  use_boot=TRUE)
  adj_q0.9 <- adjusted_surv_quantile(adj0.9, conf_int=TRUE, use_boot=TRUE)
  expect_equal(adj_q$q_surv, adj_q0.9$q_surv)
  expect_equal(adj_q$ci_lower, adj_q0.9$ci_lower)
  expect_equal(adj_q$ci_upper, adj_q0.9$ci_upper)
})

test_that("q surv, conf_int, difference, conf_level", {
  adj_q <- adjusted_surv_quantile(adj, conf_int=TRUE, conf_level=0.9,
                                  contrast="diff")
  adj_q0.9 <- adjusted_surv_quantile(adj0.9, conf_int=TRUE, contrast="diff")
  expect_equal(adj_q$diff, adj_q0.9$diff)
  expect_equal(adj_q$ci_lower, adj_q0.9$ci_lower)
  expect_equal(adj_q$ci_upper, adj_q0.9$ci_upper)
  expect_equal(adj_q$n_boot, adj_q0.9$n_boot)
})

test_that("q surv, conf_int, ratio, conf_level", {
  adj_q <- adjusted_surv_quantile(adj, conf_int=TRUE, conf_level=0.9,
                                  contrast="ratio")
  adj_q0.9 <- adjusted_surv_quantile(adj0.9, conf_int=TRUE, contrast="ratio")
  expect_equal(adj_q$ratio, adj_q0.9$ratio)
  expect_equal(adj_q$ci_lower, adj_q0.9$ci_lower)
  expect_equal(adj_q$ci_upper, adj_q0.9$ci_upper)
  expect_equal(adj_q$n_boot, adj_q0.9$n_boot)
})
