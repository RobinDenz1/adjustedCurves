library(mice)
library(survival)

set.seed(42)
sim_dat <- sim_confounded_surv(n=20)
sim_dat$group <- as.factor(sim_dat$group)

test_that("variable not in data", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                          variable="grouup",
                                                          ev_time="time",
                                                          event="event",
                                                          method="km",
                                                          conf_int=FALSE,
                                                          conf_level=0.95,
                                                          times=NULL,
                                                          bootstrap=FALSE,
                                                          n_boot=2,
                                                          na.action="na.omit",
                                                          clean_data=TRUE),
               "grouup is not a valid column name in 'data'.")
})


test_that("variable has wrong type", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                          variable="x1",
                                                          ev_time="time",
                                                          event="event",
                                                          method="km",
                                                          conf_int=FALSE,
                                                          conf_level=0.95,
                                                          times=NULL,
                                                          bootstrap=FALSE,
                                                          n_boot=2,
                                                          na.action="na.omit",
                                                          clean_data=TRUE),
              paste0("The column in 'data' specified by 'variable' needs to be",
                     " a factor variable if method='km'."))
})

test_that("variable has wrong length", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                         variable=c("x1", "x2"),
                                                         ev_time="time",
                                                         event="event",
                                                         method="km",
                                                         conf_int=FALSE,
                                                         conf_level=0.95,
                                                         times=NULL,
                                                         bootstrap=FALSE,
                                                         n_boot=2,
                                                         na.action="na.omit",
                                                         clean_data=TRUE),
                paste0("'variable' must be a character string of length 1",
                       " specifying the grouping variable in 'data'."))
})

test_that("non-standard evaluation in variable", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                          variable=group,
                                                          ev_time="time",
                                                          event="event",
                                                          method="km",
                                                          conf_int=FALSE,
                                                          conf_level=0.95,
                                                          times=NULL,
                                                          bootstrap=FALSE,
                                                          n_boot=2,
                                                          na.action="na.omit",
                                                          clean_data=TRUE),
               NULL)
})

test_that("ev_time not in data", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                          variable="group",
                                                          ev_time="time2",
                                                          event="event",
                                                          method="km",
                                                          conf_int=FALSE,
                                                          conf_level=0.95,
                                                          times=NULL,
                                                          bootstrap=FALSE,
                                                          n_boot=2,
                                                          na.action="na.omit",
                                                          clean_data=TRUE),
               "time2 is not a valid column name in 'data'.")
})

test_that("ev_time wrong format", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                          variable="group",
                                                          ev_time="group",
                                                          event="event",
                                                          method="km",
                                                          conf_int=FALSE,
                                                          conf_level=0.95,
                                                          times=NULL,
                                                          bootstrap=FALSE,
                                                          n_boot=2,
                                                          na.action="na.omit",
                                                          clean_data=TRUE),
               "The column in 'data' specified by 'ev_time' must be numeric.")
})

test_that("ev_time wrong length", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                    variable="group",
                                                    ev_time=c("group", "group"),
                                                    event="event",
                                                    method="km",
                                                    conf_int=FALSE,
                                                    conf_level=0.95,
                                                    times=NULL,
                                                    bootstrap=FALSE,
                                                    n_boot=2,
                                                    na.action="na.omit",
                                                    clean_data=TRUE),
               paste0("'ev_time' must be a character string of length 1 ",
                      "specifying the time until the event or ",
                      "censoring occured."))
})

test_that("non-standard evaluation in ev_time", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                          variable="group",
                                                          ev_time=time,
                                                          event="event",
                                                          method="km",
                                                          conf_int=FALSE,
                                                          conf_level=0.95,
                                                          times=NULL,
                                                          bootstrap=FALSE,
                                                          n_boot=2,
                                                          na.action="na.omit",
                                                          clean_data=TRUE),
               paste0("Arguments 'variable', 'ev_time', 'event' and ",
                      "'method' must be character strings, ",
                      "specifying variables in 'data'."))
})

test_that("event not in data", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                          variable="group",
                                                          ev_time="time",
                                                          event="event2",
                                                          method="km",
                                                          conf_int=FALSE,
                                                          conf_level=0.95,
                                                          times=NULL,
                                                          bootstrap=FALSE,
                                                          n_boot=2,
                                                          na.action="na.omit",
                                                          clean_data=TRUE),
               "event2 is not a valid column name in 'data'.")
})

test_that("event wrong format", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                          variable="group",
                                                          ev_time="time",
                                                          event="group",
                                                          method="km",
                                                          conf_int=FALSE,
                                                          conf_level=0.95,
                                                          times=NULL,
                                                          bootstrap=FALSE,
                                                          n_boot=2,
                                                          na.action="na.omit",
                                                          clean_data=TRUE),
      "The column in 'data' specified by 'event' must be numeric or logical.")
})

test_that("event wrong length", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                      variable="group",
                                                      ev_time="time",
                                                      event=c("group", "group"),
                                                      method="km",
                                                      conf_int=FALSE,
                                                      conf_level=0.95,
                                                      times=NULL,
                                                      bootstrap=FALSE,
                                                      n_boot=2,
                                                      na.action="na.omit",
                                                      clean_data=TRUE),
               paste0("'event' must be a character string of length 1",
                      " specifying the binary event indicator."))
})

test_that("non-standard evaluation in event", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                          variable="group",
                                                          ev_time="time",
                                                          event=event,
                                                          method="km",
                                                          conf_int=FALSE,
                                                          conf_level=0.95,
                                                          times=NULL,
                                                          bootstrap=FALSE,
                                                          n_boot=2,
                                                          na.action="na.omit",
                                                          clean_data=TRUE),
               paste0("object 'event' not found"))
})

test_that("method undefined", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                          variable="group",
                                                          ev_time="time",
                                                          event="event",
                                                          method="km3",
                                                          conf_int=FALSE,
                                                          conf_level=0.95,
                                                          times=NULL,
                                                          bootstrap=FALSE,
                                                          n_boot=2,
                                                          na.action="na.omit",
                                                          clean_data=TRUE),
               paste0("Method 'km3' is undefined. See documentation",
                      " for details on available methods."))
})

test_that("method wrong length", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                          variable="group",
                                                          ev_time="time",
                                                          event="event",
                                                          method=c("km", "km"),
                                                          conf_int=FALSE,
                                                          conf_level=0.95,
                                                          times=NULL,
                                                          bootstrap=FALSE,
                                                          n_boot=2,
                                                          na.action="na.omit",
                                                          clean_data=TRUE),
               paste0("'method' must be a single character string. ",
                      "Using multiple methods in one call is ",
                      "currently not supported."))
})

test_that("wrong conf_int type", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                          variable="group",
                                                          ev_time="time",
                                                          event="event",
                                                          method="km",
                                                          conf_int=1,
                                                          conf_level=0.95,
                                                          times=NULL,
                                                          bootstrap=FALSE,
                                                          n_boot=2,
                                                          na.action="na.omit",
                                                          clean_data=TRUE),
               "'conf_int' must be either TRUE or FALSE.")
})

test_that("wrong conf_int length", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                         variable="group",
                                                         ev_time="time",
                                                         event="event",
                                                         method="km",
                                                         conf_int=c(TRUE, TRUE),
                                                         conf_level=0.95,
                                                         times=NULL,
                                                         bootstrap=FALSE,
                                                         n_boot=2,
                                                         na.action="na.omit",
                                                         clean_data=TRUE),
               "'conf_int' must be either TRUE or FALSE, not a vector.")
})

test_that("wrong conf_level", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                          variable="group",
                                                          ev_time="time",
                                                          event="event",
                                                          method="km",
                                                          conf_int=TRUE,
                                                          conf_level=-0.95,
                                                          times=NULL,
                                                          bootstrap=FALSE,
                                                          n_boot=2,
                                                          na.action="na.omit",
                                                          clean_data=TRUE),
               "'conf_level' must be a number < 1 and > 0.")
})

test_that("wrong bootstrap type", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                          variable="group",
                                                          ev_time="time",
                                                          event="event",
                                                          method="km",
                                                          conf_int=TRUE,
                                                          conf_level=0.95,
                                                          times=NULL,
                                                          bootstrap=1,
                                                          n_boot=2,
                                                          na.action="na.omit",
                                                          clean_data=TRUE),
               "'bootstrap' must be either TRUE or FALSE.")
})

test_that("wrong bootstrap length", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                        variable="group",
                                                        ev_time="time",
                                                        event="event",
                                                        method="km",
                                                        conf_int=TRUE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=c(TRUE, TRUE),
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE),
               "'bootstrap' must be either TRUE or FALSE, not a vector.")
})

test_that("wrong n_boot type", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                          variable="group",
                                                          ev_time="time",
                                                          event="event",
                                                          method="km",
                                                          conf_int=TRUE,
                                                          conf_level=0.95,
                                                          times=NULL,
                                                          bootstrap=TRUE,
                                                          n_boot=-2,
                                                          na.action="na.omit",
                                                          clean_data=TRUE),
               "'n_boot' must be a positive integer > 2.")
})

test_that("wrong n_boot length", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                          variable="group",
                                                          ev_time="time",
                                                          event="event",
                                                          method="km",
                                                          conf_int=TRUE,
                                                          conf_level=0.95,
                                                          times=NULL,
                                                          bootstrap=TRUE,
                                                          n_boot=c(10, 10),
                                                          na.action="na.omit",
                                                          clean_data=TRUE),
               "'n_boot' must be a positive integer > 2, not a vector.")
})

test_that("wrong times", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                          variable="group",
                                                          ev_time="time",
                                                          event="event",
                                                          method="km",
                                                          conf_int=TRUE,
                                                          conf_level=0.95,
                                                          times="string",
                                                          bootstrap=TRUE,
                                                          n_boot=2,
                                                          na.action="na.omit",
                                                          clean_data=TRUE),
               "'times' must be a numeric vector or NULL.")
})

test_that("wrong na.action type", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                          variable="group",
                                                          ev_time="time",
                                                          event="event",
                                                          method="km",
                                                          conf_int=TRUE,
                                                          conf_level=0.95,
                                                          times=NULL,
                                                          bootstrap=TRUE,
                                                          n_boot=2,
                                                          na.action=1,
                                                          clean_data=TRUE),
  "'na.action' must be a function or a single character string. See details.")
})

test_that("wrong na.action length", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                          variable="group",
                                                          ev_time="time",
                                                          event="event",
                                                          method="km",
                                                          conf_int=TRUE,
                                                          conf_level=0.95,
                                                          times=NULL,
                                                          bootstrap=TRUE,
                                                          n_boot=2,
                                                          na.action=c("1", "1"),
                                                          clean_data=TRUE),
               paste0("'na.action' must be a function or a character string,",
                      " not a vector. See documentation."))
})

test_that("wrong clean_data type", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                          variable="group",
                                                          ev_time="time",
                                                          event="event",
                                                          method="km",
                                                          conf_int=TRUE,
                                                          conf_level=0.95,
                                                          times=NULL,
                                                          bootstrap=TRUE,
                                                          n_boot=2,
                                                          na.action="na.omit",
                                                          clean_data=1),
               "'clean_data' must be either TRUE or FALSE.")
})

test_that("wrong clean_data length", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                     variable="group",
                                                     ev_time="time",
                                                     event="event",
                                                     method="km",
                                                     conf_int=TRUE,
                                                     conf_level=0.95,
                                                     times=NULL,
                                                     bootstrap=TRUE,
                                                     n_boot=2,
                                                     na.action="na.omit",
                                                     clean_data=c(TRUE, TRUE)),
               "'clean_data' must be either TRUE or FALSE, not a vector.")
})

test_that("no extrapolation", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                          variable="group",
                                                          ev_time="time",
                                                          event="event",
                                                          method="km",
                                                          conf_int=TRUE,
                                                          conf_level=0.95,
                                                          times=30,
                                                          bootstrap=TRUE,
                                                          n_boot=2,
                                                          na.action="na.omit",
                                                          clean_data=TRUE),
               NULL)
})

test_that("bootstrap only with treatment_model that can be updated", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                   variable="group",
                                                   ev_time="time",
                                                   event="event",
                                                   method="iptw_km",
                                                   conf_int=FALSE,
                                                   conf_level=0.95,
                                                   times=NULL,
                                                   bootstrap=TRUE,
                                                   n_boot=2,
                                                   na.action="na.omit",
                                                   clean_data=TRUE,
                                                   treatment_model=runif(n=20)),
               paste0("'treatment_model' needs to be a model that can be ",
                      "refit or a formula object when using bootstrap=TRUE."))
})

test_that("wrong outcome_vars argument", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                        variable="group",
                                                        ev_time="time",
                                                        event="event",
                                                        method="direct_pseudo",
                                                        conf_int=FALSE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=TRUE,
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        outcome_vars=2),
               paste0("'outcome_vars' should be a character vector of ",
                      "column names in 'data', used to model ",
                      "the outcome mechanism."))
})

test_that("wrong type_time argument", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                        variable="group",
                                                        ev_time="time",
                                                        event="event",
                                                        method="direct_pseudo",
                                                        conf_int=FALSE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=TRUE,
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        outcome_vars=c("x1"),
                                                        type_time="aha"),
               "'type_time' should be either 'factor', 'bs' or 'ns'.")
})

test_that("too many spline_df", {
  expect_warning(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                         variable="group",
                                                         ev_time="time",
                                                         event="event",
                                                         method="direct_pseudo",
                                                         conf_int=FALSE,
                                                         conf_level=0.95,
                                                         times=c(0.01, 0.02),
                                                         bootstrap=TRUE,
                                                         n_boot=2,
                                                         na.action="na.omit",
                                                         clean_data=TRUE,
                                                         outcome_vars=c("x1"),
                                                         type_time="bs",
                                                         spline_df=10),
               NULL)
})

test_that("not enough time points", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                         variable="group",
                                                         ev_time="time",
                                                         event="event",
                                                         method="direct_pseudo",
                                                         conf_int=FALSE,
                                                         conf_level=0.95,
                                                         times=0.2,
                                                         bootstrap=TRUE,
                                                         n_boot=2,
                                                         na.action="na.omit",
                                                         clean_data=TRUE,
                                                         outcome_vars=c("x1")),
               paste0("'geese' models require at least two distinct ",
                      "time points. Add more points in time to 'times' ",
                      "and try again."))
})

test_that("continuous times in tmle", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                          variable="group",
                                                          ev_time="time",
                                                          event="event",
                                                          method="tmle",
                                                          conf_int=TRUE,
                                                          conf_level=0.95,
                                                          times=0.2,
                                                          bootstrap=TRUE,
                                                          n_boot=2,
                                                          na.action="na.omit",
                                                          clean_data=TRUE,
                                                          outcome_vars=c("x1")),
               paste0("Only integer time is allowed when ",
                      "using method='tmle' or method='ostmle'."))
})

test_that("wrong treatment_vars argument", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                        variable="group",
                                                        ev_time="time",
                                                        event="event",
                                                        method="emp_lik",
                                                        conf_int=TRUE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=TRUE,
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        treatment_vars=1),
               paste0("'treatment_vars' should be a character vector of ",
                      "column names in 'data', used to model the ",
                      "outcome mechanism."))
})

test_that("wrong moment argument", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                        variable="group",
                                                        ev_time="time",
                                                        event="event",
                                                        method="emp_lik",
                                                        conf_int=TRUE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=TRUE,
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        treatment_vars=c("x1"),
                                                        moment=1),
               "Argument 'moment' must be either 'first' or 'second'.")
})

test_that("0/1 variables in method='emp_lik'", {
  expect_warning(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                        variable="group",
                                                        ev_time="time",
                                                        event="event",
                                                        method="emp_lik",
                                                        conf_int=FALSE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=TRUE,
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        treatment_vars=c("x1")),
               paste0("Dichotomous variables coded with 0 and 1 found in",
                      "  'treatment_vars'. Consider recoding to -1 and 1 ",
                      "to avoid estimation problems."))
})

test_that("weights in matching", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                     variable="group",
                                     ev_time="time",
                                     event="event",
                                     method="matching",
                                     conf_int=FALSE,
                                     conf_level=0.95,
                                     times=0.2,
                                     bootstrap=FALSE,
                                     n_boot=2,
                                     na.action="na.omit",
                                     clean_data=TRUE,
                                     treatment_model=runif(20, min=10, max=20)),
                 paste0("Propensity Scores > 1 or < 0 not allowed. ",
                        "Perhaps you supplied weights on accident?"))
})

test_that("no models with method='aiptw'", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                          variable="group",
                                                          ev_time="time",
                                                          event="event",
                                                          method="aiptw",
                                                          conf_int=TRUE,
                                                          conf_level=0.95,
                                                          times=NULL,
                                                          bootstrap=FALSE,
                                                          n_boot=2,
                                                          na.action="na.omit",
                                                          clean_data=TRUE),
    "Argument 'treatment_model' must be specified when using method='aiptw'.")
})

outcome_model <- list()
class(outcome_model) <- "pecRpart"

test_that("bootstrapping with pecRpart", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                  variable="group",
                                                  ev_time="time",
                                                  event="event",
                                                  method="direct",
                                                  conf_int=TRUE,
                                                  conf_level=0.95,
                                                  times=NULL,
                                                  bootstrap=TRUE,
                                                  n_boot=2,
                                                  na.action="na.omit",
                                                  clean_data=TRUE,
                                                  outcome_model=outcome_model),
               paste0("Bootstrapping is currently not supported with ",
                      "method='direct' and an 'outcome_model' of ",
                      "class 'pecRpart'."))
})

class(outcome_model) <- "glm"

test_that("glm with censoring", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                   variable="group",
                                                   ev_time="time",
                                                   event="event",
                                                   method="direct",
                                                   conf_int=TRUE,
                                                   conf_level=0.95,
                                                   times=NULL,
                                                   bootstrap=TRUE,
                                                   n_boot=2,
                                                   na.action="na.omit",
                                                   clean_data=TRUE,
                                                   outcome_model=outcome_model),
               NULL)
})

outcome_model <- coxph(Surv(time, event) ~ 1, data=sim_dat)

test_that("coxph not including 'variable'", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                  variable="group",
                                                  ev_time="time",
                                                  event="event",
                                                  method="direct",
                                                  conf_int=TRUE,
                                                  conf_level=0.95,
                                                  times=NULL,
                                                  bootstrap=TRUE,
                                                  n_boot=2,
                                                  na.action="na.omit",
                                                  clean_data=TRUE,
                                                  outcome_model=outcome_model),
               "'variable' has to be included in the cox-regression model.")
})

sim_dat$x7 <- stats::runif(n=nrow(sim_dat))

test_that("strat_cupples with continuous confounder", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                          variable="group",
                                                          ev_time="time",
                                                          event="event",
                                                          method="strat_cupples",
                                                          conf_int=FALSE,
                                                          conf_level=0.95,
                                                          times=NULL,
                                                          bootstrap=TRUE,
                                                          n_boot=2,
                                                          na.action="na.omit",
                                                          clean_data=TRUE,
                                                          adjust_vars="x7"),
               paste0("Variables in 'adjust_vars' have to be integer, ",
                      "factor or character variables. Continuous variables ",
                      "are not allowed when using method='strat_cupples'."))
})

test_that("adjust_vars not in data", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                      variable="group",
                                                      ev_time="time",
                                                      event="event",
                                                      method="strat_amato",
                                                      conf_int=FALSE,
                                                      conf_level=0.95,
                                                      times=NULL,
                                                      bootstrap=TRUE,
                                                      n_boot=2,
                                                      na.action="na.omit",
                                                      clean_data=TRUE,
                                                      adjust_vars="x8"),
               "'adjust_vars' must specify valid column names in 'data'.")
})

test_that("invalid reference data", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                             variable="group",
                                             ev_time="time",
                                             event="event",
                                             method="strat_amato",
                                             conf_int=FALSE,
                                             conf_level=0.95,
                                             times=NULL,
                                             bootstrap=TRUE,
                                             n_boot=2,
                                             na.action="na.omit",
                                             clean_data=TRUE,
                                             adjust_vars="x1",
                                             reference=sim_dat[,c("x2", "x3")]),
               paste0("If a 'reference' data.frame is supplied, it ",
                      "needs to contain all variables listed ",
                      "in 'adjust_vars'."))
})

sim_dat$.ALL <- 1

test_that("strat_cupples with '.ALL' column", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                   variable="group",
                                                   ev_time="time",
                                                   event="event",
                                                   method="strat_cupples",
                                                   conf_int=FALSE,
                                                   conf_level=0.95,
                                                   times=NULL,
                                                   bootstrap=TRUE,
                                                   n_boot=2,
                                                   na.action="na.omit",
                                                   clean_data=TRUE,
                                                   adjust_vars="x1"),
               paste0("The column name '.ALL' cannot be used with ",
                      "method='strat_cupples'. Please rename that ",
                      "variable and rerun the function."))
})

test_that("'adjust_vars' not specified with strat_ method", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                          variable="group",
                                                          ev_time="time",
                                                          event="event",
                                                          method="strat_amato",
                                                          conf_int=FALSE,
                                                          conf_level=0.95,
                                                          times=NULL,
                                                          bootstrap=TRUE,
                                                          n_boot=2,
                                                          na.action="na.omit",
                                                          clean_data=TRUE),
"Argument 'adjust_vars' needs to be specified when using method='strat_amato'.")
})

sim_dat$.COVARS <- 1

test_that("strat_amato with '.COVARS' column", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=sim_dat,
                                                          variable="group",
                                                          ev_time="time",
                                                          event="event",
                                                          method="strat_amato",
                                                          conf_int=FALSE,
                                                          conf_level=0.95,
                                                          times=NULL,
                                                          bootstrap=TRUE,
                                                          n_boot=2,
                                                          na.action="na.omit",
                                                          clean_data=TRUE,
                                                          adjust_vars="x1"),
               paste0("The column name '.COVARS' cannot be used with ",
                      "method='strat_amato'. Please rename that variable ",
                      "and rerun the function."))
})

## multiple imputation stuff

sim_dat_na_time <- sim_dat_na_group <- sim_dat_na_event <- sim_dat_na_x1 <-
  sim_dat

sim_dat_na_time$time[3] <- NA
sim_dat_na_group$group[3] <- NA
sim_dat_na_event$event[3] <- NA
sim_dat_na_x1$x1[3] <- NA

mids_na_time <- mice(sim_dat_na_time, m=2, printFlag=FALSE)
mids_na_group <- mice(sim_dat_na_group, m=2, printFlag=FALSE)
mids_na_event <- mice(sim_dat_na_event, m=2, printFlag=FALSE)
mids_na_x1 <- mice(sim_dat_na_x1, m=2, printFlag=FALSE)

outcome_model <- list()
class(outcome_model) <- "coxph"

treatment_model <- list()
class(treatment_model) <- "glm"

test_that("wrong outcome_model with mira", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=mids_na_x1,
                                                   variable="group",
                                                   ev_time="time",
                                                   event="event",
                                                   method="direct",
                                                   conf_int=TRUE,
                                                   conf_level=0.95,
                                                   times=NULL,
                                                   bootstrap=FALSE,
                                                   n_boot=2,
                                                   na.action="na.omit",
                                                   clean_data=TRUE,
                                                   outcome_model=outcome_model),
               paste0("When using multiple imputation, mira objects need ",
                      "to be supplied to 'outcome_model' instead of ",
                      "single models. See documentation."))
})

test_that("wrong treatment_model with mira", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=mids_na_x1,
                                               variable="group",
                                               ev_time="time",
                                               event="event",
                                               method="iptw_km",
                                               conf_int=TRUE,
                                               conf_level=0.95,
                                               times=NULL,
                                               bootstrap=FALSE,
                                               n_boot=2,
                                               na.action="na.omit",
                                               clean_data=TRUE,
                                               treatment_model=treatment_model),
               paste0("When using multiple imputation, mira objects or a ",
                      "formula need to be supplied to 'treatment_model' ",
                      "instead of single models. See documentation."))
})

test_that("wrong censoring_model with mira", {
  expect_error(adjustedCurves:::check_inputs_adjustedsurv(data=mids_na_x1,
                                                 variable="group",
                                                 ev_time="time",
                                                 event="event",
                                                 method="aiptw",
                                                 conf_int=TRUE,
                                                 conf_level=0.95,
                                                 times=NULL,
                                                 bootstrap=FALSE,
                                                 n_boot=2,
                                                 na.action="na.omit",
                                                 clean_data=TRUE,
                                                 censoring_model=outcome_model),
               paste0("When using multiple imputation, mira objects need to ",
                      "be supplied to 'censoring_model' instead of single ",
                      "models. See documentation."))
})

test_that("warning with missing values in variable", {
  expect_warning(adjustedCurves:::check_inputs_adjustedsurv(data=mids_na_group,
                                                          variable="group",
                                                          ev_time="time",
                                                          event="event",
                                                          method="km",
                                                          conf_int=TRUE,
                                                          conf_level=0.95,
                                                          times=NULL,
                                                          bootstrap=FALSE,
                                                          n_boot=2,
                                                          na.action="na.omit",
                                                          clean_data=TRUE),
               paste0("Using multiple imputation with missing values in ",
                      "'variable' has not been tested yet. Use with caution."))
})

test_that("warning with missing values in ev_time", {
  expect_warning(adjustedCurves:::check_inputs_adjustedsurv(
    data=mids_na_time,
    variable="group",
    ev_time="time",
    event="event",
    method="km",
    conf_int=TRUE,
    conf_level=0.95,
    times=NULL,
    bootstrap=FALSE,
    n_boot=2,
    na.action="na.omit",
    clean_data=TRUE),
  paste0("Using multiple imputation with missing values in 'ev_time' ",
         "variable has not been tested yet. Use with caution."))
})

test_that("warning with missing values in event", {
  expect_warning(adjustedCurves:::check_inputs_adjustedsurv(data=mids_na_event,
                                                            variable="group",
                                                            ev_time="time",
                                                            event="event",
                                                            method="km",
                                                            conf_int=TRUE,
                                                            conf_level=0.95,
                                                            times=NULL,
                                                            bootstrap=FALSE,
                                                            n_boot=2,
                                                            na.action="na.omit",
                                                            clean_data=TRUE),
                 paste0("Using multiple imputation with missing values in ",
                        "'event' variable has not been tested yet. ",
                        "Use with caution."))
})

