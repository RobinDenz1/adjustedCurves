
set.seed(42)

test_that("sim, default parameters", {
  sim_dat <- sim_confounded_surv(n=15)
  expect_true(nrow(sim_dat)==15)
  expect_true(ncol(sim_dat)==9)
  expect_true(all(apply(sim_dat, 2, is.numeric)))
  expect_true(length(unique(sim_dat$event))==2)
  expect_true(length(unique(sim_dat$group))==2)
})

test_that("sim, custom values", {
  lcovars <- list(x1=c("rbinom", 1, 0.2),
                  x2=c("rbinom", 1, 0.8),
                  x3=c("rbinom", 1, 0.3),
                  x4=c("rnorm", 0, 5, -2, 2),
                  x5=c("rnorm", 2, 1, -2, 2))
  outcome_betas <- c(x1=0.2, x2=-2, x3=0, x4=0, x5=1.2)
  treatment_betas <- c(x1=0.2, x2=-2, x3=0, x4=0, x5=1.2)
  sim_dat <- sim_confounded_surv(n=20,
                                 lcovars=lcovars,
                                 outcome_betas=outcome_betas,
                                 treatment_betas=treatment_betas,
                                 group_beta=0.2)

  expect_true(nrow(sim_dat)==20)
  expect_true(ncol(sim_dat)==8)
  expect_true(all(apply(sim_dat, 2, is.numeric)))
  expect_true(length(unique(sim_dat$event))==2)
  expect_true(length(unique(sim_dat$group))==2)
})
