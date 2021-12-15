
custom_cens <- function(n){stats::rweibull(n, 1, 2)}

lcovars <- list(x1=c("rnorm", 1, 2),
                x2=c("rnorm", 3, 4),
                x3=c("runif", 1, 2))
treatment_betas <- c(x1=0.2, x2=0.6, x3=-0.9)
outcome_betas <- c(x1=1.1, x2=0, x3=-0.3)

test_that("wrong n format", {
  expect_error(adjustedCurves:::check_inputs_sim_fun(n="aha",
                                                     lcovars=NULL,
                                                     outcome_betas=NULL,
                                                     group_beta=-1,
                                                     surv_dist="weibull",
                                                     gamma=1.8,
                                                     lambda=2,
                                                     treatment_betas=NULL,
                                                     intercept=-0.5,
                                                     gtol=0.001,
                                                     cens_fun=custom_cens,
                                                     cens_args=list(),
                                                     max_t=Inf),
               "'n' must be a single positive integer.")
})

test_that("wrong n format, double", {
  expect_error(adjustedCurves:::check_inputs_sim_fun(n=1.1,
                                                     lcovars=NULL,
                                                     outcome_betas=NULL,
                                                     group_beta=-1,
                                                     surv_dist="weibull",
                                                     gamma=1.8,
                                                     lambda=2,
                                                     treatment_betas=NULL,
                                                     intercept=-0.5,
                                                     gtol=0.001,
                                                     cens_fun=custom_cens,
                                                     cens_args=list(),
                                                     max_t=Inf),
               "'n' must be a single positive integer, not a double.")
})

test_that("wrong surv_dist", {
  expect_error(adjustedCurves:::check_inputs_sim_fun(n=10,
                                                     lcovars=NULL,
                                                     outcome_betas=NULL,
                                                     group_beta=-1,
                                                     surv_dist="aha",
                                                     gamma=1.8,
                                                     lambda=2,
                                                     treatment_betas=NULL,
                                                     intercept=-0.5,
                                                     gtol=0.001,
                                                     cens_fun=custom_cens,
                                                     cens_args=list(),
                                                     max_t=Inf),
               "'surv_dist' must be either 'weibull' or 'exponential'.")
})

test_that("wrong surv_dist 2", {
  expect_error(adjustedCurves:::check_inputs_sim_fun(n=10,
                                                     lcovars=NULL,
                                                     outcome_betas=NULL,
                                                     group_beta=-1,
                                                     surv_dist=1,
                                                     gamma=1.8,
                                                     lambda=2,
                                                     treatment_betas=NULL,
                                                     intercept=-0.5,
                                                     gtol=0.001,
                                                     cens_fun=custom_cens,
                                                     cens_args=list(),
                                                     max_t=Inf),
               "'surv_dist' must be either 'weibull' or 'exponential'.")
})

test_that("wrong gamma format", {
  expect_error(adjustedCurves:::check_inputs_sim_fun(n=10,
                                                     lcovars=NULL,
                                                     outcome_betas=NULL,
                                                     group_beta=-1,
                                                     surv_dist="weibull",
                                                     gamma="1.8",
                                                     lambda=2,
                                                     treatment_betas=NULL,
                                                     intercept=-0.5,
                                                     gtol=0.001,
                                                     cens_fun=custom_cens,
                                                     cens_args=list(),
                                                     max_t=Inf),
               "'gamma' must be a single number. See details.")
})

test_that("wrong lambda format", {
  expect_error(adjustedCurves:::check_inputs_sim_fun(n=10,
                                                     lcovars=NULL,
                                                     outcome_betas=NULL,
                                                     group_beta=-1,
                                                     surv_dist="weibull",
                                                     gamma=1.8,
                                                     lambda="2",
                                                     treatment_betas=NULL,
                                                     intercept=-0.5,
                                                     gtol=0.001,
                                                     cens_fun=custom_cens,
                                                     cens_args=list(),
                                                     max_t=Inf),
               "'lambda' must be a single number. See details.")
})

test_that("wrong intercept format", {
  expect_error(adjustedCurves:::check_inputs_sim_fun(n=10,
                                                     lcovars=NULL,
                                                     outcome_betas=NULL,
                                                     group_beta=-1,
                                                     surv_dist="weibull",
                                                     gamma=1.8,
                                                     lambda=2,
                                                     treatment_betas=NULL,
                                                     intercept="-0.5",
                                                     gtol=0.001,
                                                     cens_fun=custom_cens,
                                                     cens_args=list(),
                                                     max_t=Inf),
               "'intercept' must be a single number. See details.")
})

test_that("wrong gtol format", {
  expect_error(adjustedCurves:::check_inputs_sim_fun(n=10,
                                                     lcovars=NULL,
                                                     outcome_betas=NULL,
                                                     group_beta=-1,
                                                     surv_dist="weibull",
                                                     gamma=1.8,
                                                     lambda=2,
                                                     treatment_betas=NULL,
                                                     intercept=-0.5,
                                                     gtol="10",
                                                     cens_fun=custom_cens,
                                                     cens_args=list(),
                                                     max_t=Inf),
               "'gtol' must be a number. See details.")
})

test_that("wrong gtol format, too big", {
  expect_error(adjustedCurves:::check_inputs_sim_fun(n=10,
                                                     lcovars=NULL,
                                                     outcome_betas=NULL,
                                                     group_beta=-1,
                                                     surv_dist="weibull",
                                                     gamma=1.8,
                                                     lambda=2,
                                                     treatment_betas=NULL,
                                                     intercept=-0.5,
                                                     gtol=10,
                                                     cens_fun=custom_cens,
                                                     cens_args=list(),
                                                     max_t=Inf),
               "'gtol' must be <= 1 and >= 0. See details.")
})

test_that("wrong cens_fun format", {
  expect_error(adjustedCurves:::check_inputs_sim_fun(n=10,
                                                     lcovars=NULL,
                                                     outcome_betas=NULL,
                                                     group_beta=-1,
                                                     surv_dist="weibull",
                                                     gamma=1.8,
                                                     lambda=2,
                                                     treatment_betas=NULL,
                                                     intercept=-0.5,
                                                     gtol=0.001,
                                                     cens_fun="runif",
                                                     cens_args=list(),
                                                     max_t=Inf),
               "'cens_fun' must be a function with the argument 'n' or NULL.")
})

test_that("wrong cens_args format", {
  expect_error(adjustedCurves:::check_inputs_sim_fun(n=10,
                                                     lcovars=NULL,
                                                     outcome_betas=NULL,
                                                     group_beta=-1,
                                                     surv_dist="weibull",
                                                     gamma=1.8,
                                                     lambda=2,
                                                     treatment_betas=NULL,
                                                     intercept=-0.5,
                                                     gtol=0.001,
                                                     cens_fun=runif,
                                                     cens_args="",
                                                     max_t=Inf),
               paste0("'cens_args' must be a named list of arguments to be ",
                      "passed to 'cens_fun' or NULL."))
})

test_that("wrong max_t length", {
  expect_error(adjustedCurves:::check_inputs_sim_fun(n=10,
                                                     lcovars=NULL,
                                                     outcome_betas=NULL,
                                                     group_beta=-1,
                                                     surv_dist="weibull",
                                                     gamma=1.8,
                                                     lambda=2,
                                                     treatment_betas=NULL,
                                                     intercept=-0.5,
                                                     gtol=0.001,
                                                     cens_fun=custom_cens,
                                                     cens_args=list(),
                                                     max_t=c(1, 1)),
               "'max_t' must be a single number.")
})

test_that("wrong max_t format", {
  expect_error(adjustedCurves:::check_inputs_sim_fun(n=10,
                                                     lcovars=NULL,
                                                     outcome_betas=NULL,
                                                     group_beta=-1,
                                                     surv_dist="weibull",
                                                     gamma=1.8,
                                                     lambda=2,
                                                     treatment_betas=NULL,
                                                     intercept=-0.5,
                                                     gtol=0.001,
                                                     cens_fun=custom_cens,
                                                     cens_args=list(),
                                                     max_t="aha"),
               "'max_t' must be a single number.")
})

test_that("only lcovars specified", {
  expect_error(adjustedCurves:::check_inputs_sim_fun(n=10,
                                                     lcovars=list(),
                                                     outcome_betas=NULL,
                                                     group_beta=-1,
                                                     surv_dist="weibull",
                                                     gamma=1.8,
                                                     lambda=2,
                                                     treatment_betas=NULL,
                                                     intercept=-0.5,
                                                     gtol=0.001,
                                                     cens_fun=custom_cens,
                                                     cens_args=list(),
                                                     max_t=Inf),
               paste0("'lcovars', 'outcome_betas' and 'treatment_betas' ",
                      "can either all be NULL to use default values, ",
                      "or have to be specified."))
})

test_that("wrong group_beta format", {
  expect_error(adjustedCurves:::check_inputs_sim_fun(n=10,
                                                     lcovars=NULL,
                                                     outcome_betas=NULL,
                                                     group_beta="-1",
                                                     surv_dist="weibull",
                                                     gamma=1.8,
                                                     lambda=2,
                                                     treatment_betas=NULL,
                                                     intercept=-0.5,
                                                     gtol=0.001,
                                                     cens_fun=custom_cens,
                                                     cens_args=list(),
                                                     max_t=1),
               "'group_beta' must be a number.")
})

test_that("one of lcovars / outcome_betas / treatment_betas not specified", {
  expect_error(adjustedCurves:::check_inputs_sim_fun(n=10,
                                            lcovars=list(x1=c("rnorm", 1, 2),
                                                         x2=c("rnorm", 3, 4)),
                                            outcome_betas=NULL,
                                            group_beta=-1,
                                            surv_dist="weibull",
                                            gamma=1.8,
                                            lambda=2,
                                            treatment_betas=treatment_betas,
                                            intercept=-0.5,
                                            gtol=0.001,
                                            cens_fun=custom_cens,
                                            cens_args=list(),
                                            max_t=Inf),
               paste0("'lcovars', 'outcome_betas' and 'treatment_betas' ",
                      "can either all be NULL to use default values, or ",
                      "have to be specified."))
})

test_that("lcovars wrong length", {
  expect_error(adjustedCurves:::check_inputs_sim_fun(n=10,
                                              lcovars=list(x1=c("rnorm", 1, 2),
                                                           x2=c("rnorm", 3, 4)),
                                              outcome_betas=outcome_betas,
                                              group_beta=-1,
                                              surv_dist="weibull",
                                              gamma=1.8,
                                              lambda=2,
                                              treatment_betas=treatment_betas,
                                              intercept=-0.5,
                                              gtol=0.001,
                                              cens_fun=custom_cens,
                                              cens_args=list(),
                                              max_t=Inf),
  "'outcome_betas', 'treatment_betas' and 'lcovars' must have the same length.")
})

test_that("lcovars with undefined distribution", {
  expect_error(adjustedCurves:::check_inputs_sim_fun(n=10,
                                          lcovars=list(x1=c("rnorm", 1, 2),
                                                       x2=c("rnorm", 3, 4),
                                                       x3=c("undefined", 1, 2)),
                                          outcome_betas=outcome_betas,
                                          group_beta=-1,
                                          surv_dist="weibull",
                                          gamma=1.8,
                                          lambda=2,
                                          treatment_betas=treatment_betas,
                                          intercept=-0.5,
                                          gtol=0.001,
                                          cens_fun=custom_cens,
                                          cens_args=list(),
                                          max_t=Inf),
               paste0("The first element of every vector in 'lcovars' must ",
                      "be either 'rbinom', 'rnorm' or 'runif', not undefined."))
})

test_that("no names in lcovars", {
  expect_error(adjustedCurves:::check_inputs_sim_fun(n=10,
                                               lcovars=list(c("rnorm", 1, 2),
                                                           c("rnorm", 3, 4),
                                                           c("rnorm", 1, 1)),
                                               outcome_betas=outcome_betas,
                                               group_beta=-1,
                                               surv_dist="weibull",
                                               gamma=1.8,
                                               lambda=2,
                                               treatment_betas=treatment_betas,
                                               intercept=-0.5,
                                               gtol=0.001,
                                               cens_fun=custom_cens,
                                               cens_args=list(),
                                               max_t=Inf),
               "Elements in the 'lcovars' list must be named.")
})

test_that("no names in treatment_betas", {
  expect_error(adjustedCurves:::check_inputs_sim_fun(n=10,
                                                    lcovars=lcovars,
                                                    outcome_betas=outcome_betas,
                                                    group_beta=-1,
                                                    surv_dist="weibull",
                                                    gamma=1.8,
                                                    lambda=2,
                                                    treatment_betas=c(0.2, 0.6,
                                                                      -0.9),
                                                    intercept=-0.5,
                                                    gtol=0.001,
                                                    cens_fun=custom_cens,
                                                    cens_args=list(),
                                                    max_t=Inf),
               "Elements in the 'treatment_betas' vector must be named.")
})

test_that("no names in outcome_betas", {
  expect_error(adjustedCurves:::check_inputs_sim_fun(n=10,
                                                lcovars=lcovars,
                                                outcome_betas=c(1.1, 0, -0.3),
                                                group_beta=-1,
                                                surv_dist="weibull",
                                                gamma=1.8,
                                                lambda=2,
                                                treatment_betas=treatment_betas,
                                                intercept=-0.5,
                                                gtol=0.001,
                                                cens_fun=custom_cens,
                                                cens_args=list(),
                                                max_t=Inf),
               "Elements in the 'outcome_betas' vector must be named.")
})

test_that("names not the same", {
  expect_error(adjustedCurves:::check_inputs_sim_fun(n=10,
                                        lcovars=lcovars,
                                        outcome_betas=c(y1=1.1, y2=0, y3=-0.3),
                                        group_beta=-1,
                                        surv_dist="weibull",
                                        gamma=1.8,
                                        lambda=2,
                                        treatment_betas=treatment_betas,
                                        intercept=-0.5,
                                        gtol=0.001,
                                        cens_fun=custom_cens,
                                        cens_args=list(),
                                        max_t=Inf),
               paste0("The names of the objects in 'lcovars', ",
                      "'outcome_betas' and  'treatment_betas' ",
                      "must be the same."))
})
