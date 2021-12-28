
set.seed(42)

surv <- stats::rnorm(n=50, mean=100)
se <- stats::rnorm(n=50, mean=100)

test_that("plain, bigger ci with higher conf_level", {
  ci_0.95 <- confint_surv(surv=surv, se=se, conf_level=0.95, conf_type="plain")
  ci_0.85 <- confint_surv(surv=surv, se=se, conf_level=0.85, conf_type="plain")

  ci_0.95_width <- ci_0.95$right - ci_0.95$left
  ci_0.85_width <- ci_0.85$right - ci_0.85$left

  expect_true(all(ci_0.95_width > ci_0.85_width))
})

test_that("log, bigger ci with higher conf_level", {
  ci_0.95 <- confint_surv(surv=surv, se=se, conf_level=0.95, conf_type="log")
  ci_0.85 <- confint_surv(surv=surv, se=se, conf_level=0.85, conf_type="log")

  ci_0.95_width <- ci_0.95$right - ci_0.95$left
  ci_0.85_width <- ci_0.85$right - ci_0.85$left

  expect_true(all(ci_0.95_width > ci_0.85_width))
})
