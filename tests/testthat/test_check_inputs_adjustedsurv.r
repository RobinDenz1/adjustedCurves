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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
                 NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
               NULL)
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
  NULL)
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
                 NULL)
})

