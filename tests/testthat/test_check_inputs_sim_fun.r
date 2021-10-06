
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
                                                     gtol=10,
                                                     cens_fun=custom_cens,
                                                     cens_args=list(),
                                                     max_t=Inf),
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
})
