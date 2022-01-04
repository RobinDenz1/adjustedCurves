
dat <- data.frame(time=c(0.1, 0.2, 0.3, 0.4, 0.5,
                         0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
                  surv=c(0.95, 0.8, 0.7, 0.6, 0.5,
                         0.99, 0.85, 0.75, 0.4, 0.1, 0.05),
                  group=c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1))
dat$group <- as.factor(dat$group)

dat_ci <- dat
dat_ci$se <- 1
dat_ci$ci_lower <- 1
dat_ci$ci_upper <- 10

dat_ci_boot <- dat_ci
dat_ci_boot$n_boot <- 10
dat_ci_boot$boot_surv <- 0.5
dat_ci_boot <- dplyr::select(dat_ci_boot, c("time", "group", "boot_surv",
                                            "se", "ci_lower", "ci_upper",
                                            "n_boot", "surv"))

test_that("needs zeros, no ci", {
  new_dat <- add_rows_with_zero(dat)
  expect_true(nrow(dat)==(nrow(new_dat)-2))
  expect_true(!0 %in% dat$time & 0 %in% new_dat$time)
})

test_that("needs zeros, with ci", {
  new_dat <- add_rows_with_zero(dat_ci)
  expect_true(nrow(dat_ci)==(nrow(new_dat)-2))
  expect_true(!0 %in% dat_ci$time & 0 %in% new_dat$time)
})

test_that("needs zeros, with ci, with boot", {
  new_dat <- add_rows_with_zero(dat_ci_boot)
  expect_true(nrow(dat_ci_boot)==(nrow(new_dat)-2))
  expect_true(!0 %in% dat_ci_boot$time & 0 %in% new_dat$time)
})
