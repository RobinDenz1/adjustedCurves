
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
  expect_equal(round(adj_diff$ci_lower[1], 3), 0.246)
  expect_equal(round(adj_diff$ci_upper[1], 3), 0.733)
  expect_equal(round(adj_diff$p_value[1], 3), 0.001)
})

test_that("cif, times", {
  adj_diff <- adjusted_curve_ratio(adjcif, times=0.8)
  expect_equal(nrow(adj_diff), 1)
  expect_equal(round(adj_diff$ratio[1], 4), 1.4194)
})

test_that("cif, times, conf_int", {
  adj_diff <- adjusted_curve_ratio(adjcif, times=0.8, conf_int=TRUE,
                                   use_boot=TRUE, paired_boot=FALSE)
  expect_equal(nrow(adj_diff), 1)
  expect_equal(round(adj_diff$ratio[1], 4), 1.4194)
  expect_equal(round(adj_diff$ci_lower[1], 4), 1.3114)
  expect_equal(round(adj_diff$ci_upper[1], 4), 1.5373)
  expect_equal(round(adj_diff$p_value[1], 4), 0)
})

# both curves come from one outcome model, so they are correlated
set.seed(3)
dat_direct <- sim_confounded_surv(n=100, max_t=1.5)
dat_direct$group <- factor(dat_direct$group)
mod_direct <- survival::coxph(survival::Surv(time, event) ~ group + x1 + x2,
                              data=dat_direct, x=TRUE)
adj_direct <- adjustedsurv(data=dat_direct,
                           variable="group",
                           ev_time="time",
                           event="event",
                           method="direct",
                           outcome_model=mod_direct,
                           times=c(0.4, 0.8),
                           bootstrap=TRUE,
                           n_boot=20)

test_that("surv, use_boot, Fieller interval with paired bootstrap", {
  adj_ratio <- adjusted_curve_ratio(adj_direct, times=c(0.4, 0.8),
                                    conf_int=TRUE, use_boot=TRUE)
  z <- stats::qnorm(0.975)
  boot_data <- adj_direct$boot_data
  boot_adj <- adj_direct$boot_adj

  for (i in 1:2) {
    t <- c(0.4, 0.8)[i]
    a_boot <- boot_data$surv[boot_data$group=="0" & boot_data$time==t]
    b_boot <- boot_data$surv[boot_data$group=="1" & boot_data$time==t]
    a <- boot_adj$surv[boot_adj$group=="0" & boot_adj$time==t]
    b <- boot_adj$surv[boot_adj$group=="1" & boot_adj$time==t]
    a_var <- stats::var(a_boot)
    b_var <- stats::var(b_boot)
    ab_cov <- stats::cov(a_boot, b_boot)

    # the Fieller limits are the roots of
    # (a - r*b)^2 = z^2 * (Var(a) - 2*r*Cov(a, b) + r^2*Var(b))
    fieller <- function(r) {
      (a - r*b)^2 - z^2 * (a_var - 2*r*ab_cov + r^2*b_var)
    }
    expect_equal(adj_ratio$ratio[i], a / b)
    expect_equal(fieller(adj_ratio$ci_lower[i]), 0, tolerance=1e-10)
    expect_equal(fieller(adj_ratio$ci_upper[i]), 0, tolerance=1e-10)
    expect_true(adj_ratio$ci_lower[i] < a / b && a / b < adj_ratio$ci_upper[i])
    expect_equal(adj_ratio$p_value[i],
                 2 * stats::pnorm(-abs((a - b) /
                                         sqrt(a_var + b_var - 2*ab_cov))))
  }
})

test_that("fieller_ratio_ci() with zero covariance is unchanged", {
  ci_0 <- fieller_ratio_ci(a=0.6, b=0.4, a_se=0.05, b_se=0.04)
  ci_cov <- fieller_ratio_ci(a=0.6, b=0.4, a_se=0.05, b_se=0.04, ab_cov=0)
  expect_identical(ci_0, ci_cov)
  expect_equal(round(ci_0$ci_lower, 4), 1.1651)
  expect_equal(round(ci_0$ci_upper, 4), 1.9548)

  # zero variance in b, e.g. survival probabilities at time 0
  ci_1 <- fieller_ratio_ci(a=1, b=1, a_se=0, b_se=0, ab_cov=0)
  expect_equal(c(ci_1$ci_lower, ci_1$ci_upper), c(1, 1))
})
