
custom_cens <- function(n) {
  stats::rweibull(n, 1, 2)
}

outcome_betas <- list(c(0.03, 0.4),
                      c(1.1, 0.8),
                      c(0, 0))

treatment_betas <- c(x1=0, x2=log(3), x3=log(1.2))

lcovars <- list(x1=c("rbinom", 1, 0.3),
                x2=c("rbinom", 1, 0.7),
                x3=c("rbinom", 1, 0.5))

test_that("wrong n format", {
  expect_error(check_inputs_sim_crisk_fun(n="aha",
                                          lcovars=NULL,
                                          outcome_betas=NULL,
                                          group_beta=c(1, 0),
                                          gamma=c(1.8, 1.8),
                                          lambda=c(2, 2),
                                          treatment_betas=NULL,
                                          intercept=-0.5,
                                          gtol=0.001,
                                          cens_fun=custom_cens,
                                          cens_args=list(),
                                          max_t=1.7),
               "'n' must be a single positive integer.")
})

test_that("wrong n format double", {
  expect_error(check_inputs_sim_crisk_fun(n=1.3,
                                          lcovars=NULL,
                                          outcome_betas=NULL,
                                          group_beta=c(1, 0),
                                          gamma=c(1.8, 1.8),
                                          lambda=c(2, 2),
                                          treatment_betas=NULL,
                                          intercept=-0.5,
                                          gtol=0.001,
                                          cens_fun=custom_cens,
                                          cens_args=list(),
                                          max_t=1.7),
               "'n' must be a positive integer, not a double.")
})

test_that("wrong gamma format", {
  expect_error(check_inputs_sim_crisk_fun(n=10,
                                          lcovars=NULL,
                                          outcome_betas=NULL,
                                          group_beta=c(1, 0),
                                          gamma=c("1.8", "1.8"),
                                          lambda=c(2, 2),
                                          treatment_betas=NULL,
                                          intercept=-0.5,
                                          gtol=0.001,
                                          cens_fun=custom_cens,
                                          cens_args=list(),
                                          max_t=1.7),
               "'gamma' must be a numeric vector with length > 1. See details.")
})

test_that("wrong lambda format", {
  expect_error(check_inputs_sim_crisk_fun(n=10,
                                          lcovars=NULL,
                                          outcome_betas=NULL,
                                          group_beta=c(1, 0),
                                          gamma=c(1.8, 1.8),
                                          lambda=c("2", "2"),
                                          treatment_betas=NULL,
                                          intercept=-0.5,
                                          gtol=0.001,
                                          cens_fun=custom_cens,
                                          cens_args=list(),
                                          max_t=1.7),
              "'lambda' must be a numeric vector with length > 1. See details.")
})

test_that("wrong intercept format", {
  expect_error(check_inputs_sim_crisk_fun(n=10,
                                          lcovars=NULL,
                                          outcome_betas=NULL,
                                          group_beta=c(1, 0),
                                          gamma=c(1.8, 1.8),
                                          lambda=c(2, 2),
                                          treatment_betas=NULL,
                                          intercept="-0.5",
                                          gtol=0.001,
                                          cens_fun=custom_cens,
                                          cens_args=list(),
                                          max_t=1.7),
               "'intercept' must be a number. See details.")
})

test_that("wrong gtol too big", {
  expect_error(check_inputs_sim_crisk_fun(n=10,
                                          lcovars=NULL,
                                          outcome_betas=NULL,
                                          group_beta=c(1, 0),
                                          gamma=c(1.8, 1.8),
                                          lambda=c(2, 2),
                                          treatment_betas=NULL,
                                          intercept=-0.5,
                                          gtol=10,
                                          cens_fun=custom_cens,
                                          cens_args=list(),
                                          max_t=1.7),
               "'gtol' must be <= 1 and >= 0. See details.")
})

test_that("wrong gtol format", {
  expect_error(check_inputs_sim_crisk_fun(n=10,
                                          lcovars=NULL,
                                          outcome_betas=NULL,
                                          group_beta=c(1, 0),
                                          gamma=c(1.8, 1.8),
                                          lambda=c(2, 2),
                                          treatment_betas=NULL,
                                          intercept=-0.5,
                                          gtol="0.1",
                                          cens_fun=custom_cens,
                                          cens_args=list(),
                                          max_t=1.7),
               "'gtol' must be a number. See details.")
})

test_that("wrong cens_fun format", {
  expect_error(check_inputs_sim_crisk_fun(n=10,
                                          lcovars=NULL,
                                          outcome_betas=NULL,
                                          group_beta=c(1, 0),
                                          gamma=c(1.8, 1.8),
                                          lambda=c(2, 2),
                                          treatment_betas=NULL,
                                          intercept=-0.5,
                                          gtol=0.001,
                                          cens_fun="runif",
                                          cens_args=list(),
                                          max_t=1.7),
              "'cens_fun' must be a function with the argument 'n' or NULL.")
})

test_that("wrong cens_args format", {
  expect_error(check_inputs_sim_crisk_fun(n=10,
                                          lcovars=NULL,
                                          outcome_betas=NULL,
                                          group_beta=c(1, 0),
                                          gamma=c(1.8, 1.8),
                                          lambda=c(2, 2),
                                          treatment_betas=NULL,
                                          intercept=-0.5,
                                          gtol=0.001,
                                          cens_fun=runif,
                                          cens_args="list()",
                                          max_t=1.7),
              paste0("'cens_args' must be a named list of arguments ",
                     "to be passed to 'cens_fun' or NULL."))
})

test_that("wrong max_t format", {
  expect_error(check_inputs_sim_crisk_fun(n=10,
                                          lcovars=NULL,
                                          outcome_betas=NULL,
                                          group_beta=c(1, 0),
                                          gamma=c(1.8, 1.8),
                                          lambda=c(2, 2),
                                          treatment_betas=NULL,
                                          intercept=-0.5,
                                          gtol=0.001,
                                          cens_fun=custom_cens,
                                          cens_args=list(),
                                          max_t="1.7"),
               "'max_t' must be a number.")
})

test_that("max_t too small", {
  expect_error(check_inputs_sim_crisk_fun(n=10,
                                          lcovars=NULL,
                                          outcome_betas=NULL,
                                          group_beta=c(1, 0),
                                          gamma=c(1.8, 1.8),
                                          lambda=c(2, 2),
                                          treatment_betas=NULL,
                                          intercept=-0.5,
                                          gtol=0.001,
                                          cens_fun=custom_cens,
                                          cens_args=list(),
                                          max_t=0),
               "'max_t' must be bigger than zero.")
})

test_that("wrong group_beta format", {
  expect_error(check_inputs_sim_crisk_fun(n=10,
                                          lcovars=NULL,
                                          outcome_betas=NULL,
                                          group_beta=c("1", "0"),
                                          gamma=c(1.8, 1.8),
                                          lambda=c(2, 2),
                                          treatment_betas=NULL,
                                          intercept=-0.5,
                                          gtol=0.001,
                                          cens_fun=custom_cens,
                                          cens_args=list(),
                                          max_t=1),
               "'group_beta' must be a numeric vector with length > 1.")
})

test_that("wrong max_t length", {
  expect_error(check_inputs_sim_crisk_fun(n=10,
                                          lcovars=NULL,
                                          outcome_betas=NULL,
                                          group_beta=c(1, 0),
                                          gamma=c(1.8, 1.8),
                                          lambda=c(2, 2),
                                          treatment_betas=NULL,
                                          intercept=-0.5,
                                          gtol=0.001,
                                          cens_fun=custom_cens,
                                          cens_args=list(),
                                          max_t=c(1, 1.2)),
               "'max_t' must be a number.")
})

test_that("only lcovars specified", {
  expect_error(check_inputs_sim_crisk_fun(n=10,
                                          lcovars=list(),
                                          outcome_betas=NULL,
                                          group_beta=c(1, 0),
                                          gamma=c(1.8, 1.8),
                                          lambda=c(2, 2),
                                          treatment_betas=NULL,
                                          intercept=-0.5,
                                          gtol=0.001,
                                          cens_fun=custom_cens,
                                          cens_args=list(),
                                          max_t=1.7),
               paste0("'lcovars', 'outcome_betas' and 'treatment_betas' ",
                      "can either all be NULL to use default values, or ",
                      "have to be specified."))
})

test_that("lcovars wrong length", {
  expect_error(check_inputs_sim_crisk_fun(n=10,
                                          lcovars=list(x1=c("rbinom", 1, 0.3),
                                                       x2=c("rbinom", 1, 0.7)),
                                          outcome_betas=outcome_betas,
                                          group_beta=c(1, 0),
                                          gamma=c(1.8, 1.8),
                                          lambda=c(2, 2),
                                          treatment_betas=treatment_betas,
                                          intercept=-0.5,
                                          gtol=0.001,
                                          cens_fun=custom_cens,
                                          cens_args=list(),
                                          max_t=1.7),
  "'outcome_betas', 'treatment_betas' and 'lcovars' must have the same length.")
})

test_that("lcovars with undefined distribution", {
  expect_error(check_inputs_sim_crisk_fun(n=10,
                                          lcovars=list(x1=c("rnorm", 1, 2),
                                                       x2=c("rnorm", 3, 4),
                                                       x3=c("undefined", 1, 2)),
                                          outcome_betas=outcome_betas,
                                          group_beta=c(1, 0),
                                          gamma=c(1.8, 1.8),
                                          lambda=c(2, 2),
                                          treatment_betas=treatment_betas,
                                          intercept=-0.5,
                                          gtol=0.001,
                                          cens_fun=custom_cens,
                                          cens_args=list(),
                                          max_t=1.7),
               paste0("The first element of every vector in 'lcovars' ",
                      "must be either 'rbinom', 'rnorm' or 'runif', ",
                      "not undefined."))
})

test_that("no names in lcovars", {
  expect_error(check_inputs_sim_crisk_fun(n=10,
                                          lcovars=list(c("rnorm", 1, 2),
                                                       c("rnorm", 3, 4),
                                                       c("rnorm", 1, 1)),
                                          outcome_betas=outcome_betas,
                                          group_beta=c(1, 0),
                                          gamma=c(1.8, 1.8),
                                          lambda=c(2, 2),
                                          treatment_betas=treatment_betas,
                                          intercept=-0.5,
                                          gtol=0.001,
                                          cens_fun=custom_cens,
                                          cens_args=list(),
                                          max_t=1.7),
               "Elements in the 'lcovars' list must be named.")
})

test_that("no names in treatment_betas", {
  expect_error(check_inputs_sim_crisk_fun(n=10,
                                         lcovars=lcovars,
                                         outcome_betas=outcome_betas,
                                         group_beta=c(1, 0),
                                         gamma=c(1.8, 1.8),
                                         lambda=c(2, 2),
                                         treatment_betas=c(0, log(3), log(1.2)),
                                         intercept=-0.5,
                                         gtol=0.001,
                                         cens_fun=custom_cens,
                                         cens_args=list(),
                                         max_t=1.7),
               "Elements in the 'treatment_betas' vector must be named.")
})

test_that("different lengths in outcome_betas", {
  expect_error(check_inputs_sim_crisk_fun(n=10,
                                          lcovars=lcovars,
                                          outcome_betas=list(c(0.03, 0.4),
                                                             c(1.1, 0.8),
                                                             c(0)),
                                          group_beta=c(1, 0),
                                          gamma=c(1.8, 1.8),
                                          lambda=c(2, 2),
                                          treatment_betas=treatment_betas,
                                          intercept=-0.5,
                                          gtol=0.001,
                                          cens_fun=custom_cens,
                                          cens_args=list(),
                                          max_t=1.7),
  "The vectors supplied in 'outcome_betas' all need to have, the same length.")
})

test_that("different lengths in group_beta / gamma / lambda", {
  expect_error(check_inputs_sim_crisk_fun(n=10,
                                          lcovars=lcovars,
                                          outcome_betas=list(c(0.03, 0.4),
                                                             c(1.1, 0.8),
                                                             c(0, 0)),
                                          group_beta=c(1, 0),
                                          gamma=c(1.8, 2, 2),
                                          lambda=c(2, 2),
                                          treatment_betas=c(x1=0,
                                                            x2=log(3),
                                                            x3=0),
                                          intercept=-0.5,
                                          gtol=0.001,
                                          cens_fun=custom_cens,
                                          cens_args=list(),
                                          max_t=1.7),
  "Arguments 'group_beta', 'gamma' and 'lambda' need to have the same length.")
})
