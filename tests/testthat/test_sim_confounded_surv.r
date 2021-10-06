library(survival)

test_that("sim, default parameters", {
  expect_error(sim_confounded_surv(n=10), NA)
})

test_that("sim, correct n", {
  expect_true(nrow(sim_confounded_surv(10))==10)
})

test_that("sim, custom values", {
  expect_error({
    lcovars <- list(x1=c("rbinom", 1, 0.2),
                    x2=c("rbinom", 1, 0.8),
                    x3=c("rbinom", 1, 0.3),
                    x4=c("rnorm", 0, 5, -2, 2),
                    x5=c("rnorm", 2, 1, -2, 2))
    outcome_betas <- c(x1=0.2, x2=-2, x3=0, x4=0, x5=1.2)
    treatment_betas <- c(x1=0.2, x2=-2, x3=0, x4=0, x5=1.2)
    sim_confounded_surv(n=100, lcovars=lcovars, outcome_betas=outcome_betas,
                        treatment_betas=treatment_betas, group_beta=0.2)
  }, NA)
})
