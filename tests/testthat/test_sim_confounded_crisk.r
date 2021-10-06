library(survival)

test_that("sim, default parameters", {
  expect_error(sim_confounded_crisk(n=10), NA)
})

test_that("sim, correct n", {
  expect_true(nrow(sim_confounded_crisk(10))==10)
})

test_that("sim, custom values", {
  expect_error({
    outcome_betas <- list(c(0.03, 0.4),
                          c(1.1, 0.8),
                          c(0, 0),
                          c(-0.2, -0.4),
                          c(log(1.3), log(1.3)/3),
                          c(0, 0))

    treatment_betas <- c(x1=0, x2=log(3), x3=log(1.2),
                         x4=0, x5=log(1.1), x6=log(1.4))

    lcovars <- list(x1=c("rbinom", 1, 0.3),
                    x2=c("rbinom", 1, 0.7),
                    x3=c("rbinom", 1, 0.5),
                    x4=c("rnorm", 0, 1),
                    x5=c("rnorm", 0, 1.1),
                    x6=c("rnorm", 0, 0.9))

    sim_dat <- sim_confounded_crisk(n=100,
                                    treatment_betas=treatment_betas,
                                    outcome_betas=outcome_betas,
                                    lcovars=lcovars)
  }, NA)
})
