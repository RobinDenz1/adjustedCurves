
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
  adj_diff <- adjusted_curve_diff(adj)
  expect_equal(nrow(adj_diff), 193)
  expect_equal(round(adj_diff$diff[1], 4), 0)
  expect_equal(round(adj_diff$diff[139], 4), -0.291)
})

test_that("surv, times", {
  adj_diff <- adjusted_curve_diff(adj, times=0.8)
  expect_equal(nrow(adj_diff), 1)
  expect_equal(round(adj_diff$diff[1], 4), -0.2425)
})

test_that("surv, times, conf_int", {
  adj_diff <- adjusted_curve_diff(adj, times=0.8, conf_int=TRUE)
  expect_equal(nrow(adj_diff), 1)
  expect_equal(round(adj_diff$diff[1], 4), -0.2425)
  expect_equal(round(adj_diff$se[1], 4), 0.0752)
  expect_equal(round(adj_diff$ci_lower[1], 3), -0.390)
})

test_that("cif, times", {
  adj_diff <- adjusted_curve_diff(adjcif, times=0.8)
  expect_equal(nrow(adj_diff), 1)
  expect_equal(round(adj_diff$diff[1], 4), 0.2329)
})

test_that("cif, times, conf_int", {
  adj_diff <- adjusted_curve_diff(adjcif, times=0.8, conf_int=TRUE,
                                  use_boot=TRUE, paired_boot=FALSE)
  expect_equal(nrow(adj_diff), 1)
  expect_equal(round(adj_diff$diff[1], 4), 0.2329)
  expect_equal(round(adj_diff$se[1], 4), 0.0271)
})

# SE of the paired bootstrap differences at time t, using the last
# estimate at or before t in each bootstrap replication
paired_boot_se <- function(boot_data, est, group_1, group_2, t) {
  boot_data <- boot_data[boot_data$time <= t & !is.na(boot_data[, est]), ]
  boot_data <- boot_data[order(boot_data$boot, boot_data$time), ]
  last_est <- function(group) {
    dat <- boot_data[boot_data$group==group, ]
    dat <- dat[!duplicated(dat$boot, fromLast=TRUE), ]
    stats::setNames(dat[, est], dat$boot)
  }
  est_1 <- last_est(group_1)
  est_2 <- last_est(group_2)
  ids <- intersect(names(est_1), names(est_2))
  stats::sd(est_1[ids] - est_2[ids])
}

test_that("cif, times, conf_int, paired bootstrap", {
  adj_diff <- adjusted_curve_diff(adjcif, times=0.8, conf_int=TRUE,
                                  use_boot=TRUE)
  expect_equal(round(adj_diff$diff[1], 4), 0.2329)
  expect_equal(adj_diff$se, paired_boot_se(adjcif$boot_data, est="cif",
                                           group_1="0", group_2="1", t=0.8))
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

test_that("surv, use_boot, paired bootstrap", {
  adj_diff <- adjusted_curve_diff(adj_direct, times=c(0.4, 0.8),
                                  conf_int=TRUE, use_boot=TRUE)
  se_paired <- c(paired_boot_se(adj_direct$boot_data, est="surv",
                                group_1="0", group_2="1", t=0.4),
                 paired_boot_se(adj_direct$boot_data, est="surv",
                                group_1="0", group_2="1", t=0.8))
  expect_equal(adj_diff$se, se_paired)
  expect_equal(adj_diff$ci_lower,
               adj_diff$diff - stats::qnorm(0.975) * adj_diff$se)
  expect_equal(adj_diff$p_value,
               2 * stats::pnorm(-abs(adj_diff$diff / adj_diff$se)))

  # paired_boot=FALSE pools the bootstrap SEs of the two curves
  adj_diff_pooled <- adjusted_curve_diff(adj_direct, times=c(0.4, 0.8),
                                         conf_int=TRUE, use_boot=TRUE,
                                         paired_boot=FALSE)
  boot_adj <- adj_direct$boot_adj
  expect_equal(adj_diff_pooled$diff, adj_diff$diff)
  expect_equal(adj_diff_pooled$se,
               sqrt(boot_adj$se[boot_adj$group=="0"]^2 +
                    boot_adj$se[boot_adj$group=="1"]^2))
})

test_that("surv, use_boot, paired bootstrap, between time points", {
  adj_diff <- adjusted_curve_diff(adj_direct, times=c(0.5, 0.7),
                                  conf_int=TRUE, use_boot=TRUE)
  expect_equal(adj_diff$se,
               c(paired_boot_se(adj_direct$boot_data, est="surv",
                                group_1="0", group_2="1", t=0.5),
                 paired_boot_se(adj_direct$boot_data, est="surv",
                                group_1="0", group_2="1", t=0.7)))
})

test_that("surv, use_boot, paired bootstrap, multiple imputation", {
  skip_if_not_installed("mice")

  set.seed(4)
  dat_mi <- sim_confounded_surv(n=80, max_t=1.5)
  dat_mi$group <- factor(dat_mi$group)
  dat_mi$x1 <- ifelse(stats::runif(n=80) <= 0.7, dat_mi$x1, NA)
  imp <- suppressWarnings(mice::mice(dat_mi, m=3, method="pmm",
                                     printFlag=FALSE))
  mods <- with(imp, survival::coxph(survival::Surv(time, event) ~ group +
                                      x1 + x2, x=TRUE))
  adj_mi <- adjustedsurv(data=imp,
                         variable="group",
                         ev_time="time",
                         event="event",
                         method="direct",
                         outcome_model=mods,
                         times=0.8,
                         bootstrap=TRUE,
                         n_boot=10)
  adj_diff <- adjusted_curve_diff(adj_mi, times=0.8, conf_int=TRUE,
                                  use_boot=TRUE)

  # Rubin's rule for the difference: mean within-imputation variance of
  # the paired differences plus (1 + 1/m) times the between variance
  var_w <- vapply(adj_mi$mids_analyses, function(x) {
    paired_boot_se(x$boot_data, est="surv", group_1="0", group_2="1",
                   t=0.8)^2
  }, numeric(1))
  diff_m <- vapply(adj_mi$mids_analyses, function(x) {
    x$boot_adj$surv[x$boot_adj$group=="0"] -
      x$boot_adj$surv[x$boot_adj$group=="1"]
  }, numeric(1))
  m <- length(adj_mi$mids_analyses)
  expect_equal(adj_diff$diff, mean(diff_m))
  expect_equal(adj_diff$se,
               sqrt(mean(var_w) + stats::var(diff_m) * (1 + 1/m)))
})
