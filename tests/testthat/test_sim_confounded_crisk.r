
set.seed(42)

test_that("sim, default parameters", {
  sim_dat <- sim_confounded_crisk(n=15)
  expect_true(nrow(sim_dat)==15)
  expect_true(ncol(sim_dat)==9)
  expect_true(all(apply(sim_dat, 2, is.numeric)))
  expect_true(length(unique(sim_dat$event))==3)
  expect_true(length(unique(sim_dat$group))==2)
})

test_that("sim, custom values", {
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

  sim_dat <- sim_confounded_crisk(n=15,
                                  treatment_betas=treatment_betas,
                                  outcome_betas=outcome_betas,
                                  lcovars=lcovars)
  expect_true(nrow(sim_dat)==15)
  expect_true(ncol(sim_dat)==9)
  expect_true(all(apply(sim_dat, 2, is.numeric)))
  expect_true(length(unique(sim_dat$event))==3)
  expect_true(length(unique(sim_dat$group))==2)
})
