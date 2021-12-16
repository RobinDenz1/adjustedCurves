library(mice)

set.seed(42)
sim_dat <- sim_confounded_crisk(n=20)
sim_dat$group <- as.factor(sim_dat$group)

test_that("variable not in data", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                        variable="grouup",
                                                        ev_time="time",
                                                        event="event",
                                                        method="aalen_johansen",
                                                        conf_int=FALSE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=FALSE,
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        cause=1),
               NULL)
})


test_that("variable has wrong type", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                        variable="x1",
                                                        ev_time="time",
                                                        event="event",
                                                        method="aalen_johansen",
                                                        conf_int=FALSE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=FALSE,
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        cause=1),
               NULL)
})

test_that("variable has wrong length", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                        variable=c("x1", "x2"),
                                                        ev_time="time",
                                                        event="event",
                                                        method="aalen_johansen",
                                                        conf_int=FALSE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=FALSE,
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        cause=1),
               NULL)
})

test_that("non-standard evaluation in variable", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                        variable=group,
                                                        ev_time="time",
                                                        event="event",
                                                        method="aalen_johansen",
                                                        conf_int=FALSE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=FALSE,
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        cause=1),
               NULL)
})

test_that("ev_time not in data", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                        variable="group",
                                                        ev_time="time2",
                                                        event="event",
                                                        method="aalen_johansen",
                                                        conf_int=FALSE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=FALSE,
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        cause=1),
               NULL)
})

test_that("ev_time wrong format", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                        variable="group",
                                                        ev_time="group",
                                                        event="event",
                                                        method="aalen_johansen",
                                                        conf_int=FALSE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=FALSE,
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        cause=1),
               NULL)
})

sim_dat$time2 <- sim_dat$time
sim_dat$time2[1] <- -1

test_that("ev_time includes negative values", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                        variable="group",
                                                        ev_time="time2",
                                                        event="event",
                                                        method="aalen_johansen",
                                                        conf_int=FALSE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=FALSE,
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        cause=1),
               "All values in the 'ev_time' variable must be >= 0.")
})

test_that("ev_time wrong length", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                      variable="group",
                                                      ev_time=c("time", "time"),
                                                      event="event",
                                                      method="aalen_johansen",
                                                      conf_int=FALSE,
                                                      conf_level=0.95,
                                                      times=NULL,
                                                      bootstrap=FALSE,
                                                      n_boot=2,
                                                      na.action="na.omit",
                                                      clean_data=TRUE,
                                                      cause=1),
               NULL)
})

test_that("non-standard evaluation in ev_time", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                        variable="group",
                                                        ev_time=time,
                                                        event="event",
                                                        method="aalen_johansen",
                                                        conf_int=FALSE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=FALSE,
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        cause=1),
               NULL)
})

test_that("event not in data", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                        variable="group",
                                                        ev_time="time",
                                                        event="event2",
                                                        method="aalen_johansen",
                                                        conf_int=FALSE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=FALSE,
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        cause=1),
               NULL)
})

test_that("event wrong format", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                        variable="group",
                                                        ev_time="time",
                                                        event="group",
                                                        method="aalen_johansen",
                                                        conf_int=FALSE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=FALSE,
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        cause=1),
               NULL)
})

sim_dat$event2 <- sample(c(0, 1), size=nrow(sim_dat), replace=TRUE)

test_that("event only includes two values", {
  expect_warning(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                        variable="group",
                                                        ev_time="time",
                                                        event="event2",
                                                        method="aalen_johansen",
                                                        conf_int=FALSE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=FALSE,
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        cause=1),
               paste0("It is recommended to use the 'adjustedsurv' function",
                      " when the 'event' variable is binary."))
})

test_that("event wrong length", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                      variable="group",
                                                      ev_time="time",
                                                      event=c("event", "event"),
                                                      method="aalen_johansen",
                                                      conf_int=FALSE,
                                                      conf_level=0.95,
                                                      times=NULL,
                                                      bootstrap=FALSE,
                                                      n_boot=2,
                                                      na.action="na.omit",
                                                      clean_data=TRUE,
                                                      cause=1),
               NULL)
})

test_that("non-standard evaluation in event", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                        variable="group",
                                                        ev_time="time",
                                                        event=event,
                                                        method="aalen_johansen",
                                                        conf_int=FALSE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=FALSE,
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        cause=1),
               NULL)
})

test_that("cause has wrong type", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                        variable="group",
                                                        ev_time="time",
                                                        event="event",
                                                        method="aalen_johansen",
                                                        conf_int=FALSE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=FALSE,
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        cause="1"),
               NULL)
})

test_that("more than one cause supplied", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                        variable="group",
                                                        ev_time="time",
                                                        event="event",
                                                        method="aalen_johansen",
                                                        conf_int=FALSE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=FALSE,
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        cause=c(1, 2)),
               NULL)
})

test_that("method undefined", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
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
                                                         clean_data=TRUE,
                                                         cause=1),
               NULL)
})

test_that("wrong method length", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
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
                                                         clean_data=TRUE,
                                                         cause=1),
               NULL)
})

test_that("wrong conf_int type", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                        variable="group",
                                                        ev_time="time",
                                                        event="event",
                                                        method="aalen_johansen",
                                                        conf_int=1,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=FALSE,
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        cause=1),
               NULL)
})

test_that("wrong conf_int length", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                        variable="group",
                                                        ev_time="time",
                                                        event="event",
                                                        method="aalen_johansen",
                                                        conf_int=c(TRUE, TRUE),
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=FALSE,
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        cause=1),
               NULL)
})

test_that("wrong conf_level type", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                        variable="group",
                                                        ev_time="time",
                                                        event="event",
                                                        method="aalen_johansen",
                                                        conf_int=TRUE,
                                                        conf_level=-0.95,
                                                        times=NULL,
                                                        bootstrap=FALSE,
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        cause=1),
               NULL)
})

test_that("wrong conf_level length", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                      variable="group",
                                                      ev_time="time",
                                                      event="event",
                                                      method="aalen_johansen",
                                                      conf_int=TRUE,
                                                      conf_level=c(0.95, 0.95),
                                                      times=NULL,
                                                      bootstrap=FALSE,
                                                      n_boot=2,
                                                      na.action="na.omit",
                                                      clean_data=TRUE,
                                                      cause=1),
               NULL)
})

test_that("wrong bootstrap type", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                        variable="group",
                                                        ev_time="time",
                                                        event="event",
                                                        method="aalen_johansen",
                                                        conf_int=TRUE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=1,
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        cause=1),
               NULL)
})

test_that("wrong bootstrap length", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                        variable="group",
                                                        ev_time="time",
                                                        event="event",
                                                        method="aalen_johansen",
                                                        conf_int=TRUE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=c(TRUE, TRUE),
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        cause=1),
               NULL)
})

test_that("wrong n_boot type", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                        variable="group",
                                                        ev_time="time",
                                                        event="event",
                                                        method="aalen_johansen",
                                                        conf_int=TRUE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=TRUE,
                                                        n_boot=-2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        cause=1),
               NULL)
})

test_that("wrong n_boot length", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                        variable="group",
                                                        ev_time="time",
                                                        event="event",
                                                        method="aalen_johansen",
                                                        conf_int=TRUE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=TRUE,
                                                        n_boot=c(2, 2),
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        cause=1),
               NULL)
})

test_that("wrong times", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                        variable="group",
                                                        ev_time="time",
                                                        event="event",
                                                        method="aalen_johansen",
                                                        conf_int=TRUE,
                                                        conf_level=0.95,
                                                        times="string",
                                                        bootstrap=TRUE,
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        cause=1),
               NULL)
})

test_that("wrong na.action type", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                        variable="group",
                                                        ev_time="time",
                                                        event="event",
                                                        method="aalen_johansen",
                                                        conf_int=TRUE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=TRUE,
                                                        n_boot=10,
                                                        na.action=1,
                                                        clean_data=TRUE,
                                                        cause=1),
               NULL)
})

test_that("wrong na.action length", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                              variable="group",
                                              ev_time="time",
                                              event="event",
                                              method="aalen_johansen",
                                              conf_int=TRUE,
                                              conf_level=0.95,
                                              times=NULL,
                                              bootstrap=TRUE,
                                              n_boot=10,
                                              na.action=c("na.omit", "na.omit"),
                                              clean_data=TRUE,
                                              cause=1),
               NULL)
})

test_that("wrong clean_data type", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                        variable="group",
                                                        ev_time="time",
                                                        event="event",
                                                        method="aalen_johansen",
                                                        conf_int=TRUE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=TRUE,
                                                        n_boot=10,
                                                        na.action="na.omit",
                                                        clean_data=1,
                                                        cause=1),
               NULL)
})

test_that("wrong clean_data length", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                      variable="group",
                                                      ev_time="time",
                                                      event="event",
                                                      method="aalen_johansen",
                                                      conf_int=TRUE,
                                                      conf_level=0.95,
                                                      times=NULL,
                                                      bootstrap=TRUE,
                                                      n_boot=10,
                                                      na.action="na.omit",
                                                      clean_data=c(TRUE, TRUE),
                                                      cause=1),
               NULL)
})

test_that("no extrapolation", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                        variable="group",
                                                        ev_time="time",
                                                        event="event",
                                                        method="aalen_johansen",
                                                        conf_int=TRUE,
                                                        conf_level=0.95,
                                                        times=30,
                                                        bootstrap=TRUE,
                                                        n_boot=2,
                                                        na.action="na.omit",
                                                        clean_data=TRUE,
                                                        cause=1),
               NULL)
})

test_that("bootstrap only with treatment_model that can be updated", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                    variable="group",
                                                    ev_time="time",
                                                    event="event",
                                                    method="iptw_pseudo",
                                                    conf_int=FALSE,
                                                    conf_level=0.95,
                                                    times=NULL,
                                                    bootstrap=TRUE,
                                                    n_boot=2,
                                                    na.action="na.omit",
                                                    clean_data=TRUE,
                                                    treatment_model=runif(n=20),
                                                    cause=1),
               NULL)
})

test_that("wrong outcome_vars argument", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
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
                                                         outcome_vars=2,
                                                         cause=1),
               NULL)
})

test_that("wrong type_time argument", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
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
                                                         type_time="aha",
                                                         cause=1),
               NULL)
})

test_that("too many spline_df", {
  expect_warning(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
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
                                                         spline_df=10,
                                                         cause=1),
                 NULL)
})

test_that("not enough time points", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                         variable="group",
                                                         ev_time="time",
                                                         event="event",
                                                         method="direct_pseudo",
                                                         conf_int=TRUE,
                                                         conf_level=0.95,
                                                         times=0.2,
                                                         bootstrap=TRUE,
                                                         n_boot=2,
                                                         na.action="na.omit",
                                                         clean_data=TRUE,
                                                         outcome_vars=c("x1"),
                                                         cause=1),
               NULL)
})

test_that("continuous times in tmle", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
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
                                                         outcome_vars=c("x1"),
                                                         cause=1),
               NULL)
})

test_that("weights in matching", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
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
                                      treatment_model=runif(20, min=10, max=20),
                                      cause=1),
               NULL)
})

test_that("no treatment_model in matching", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
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
                                      cause=1),
               NULL)
})

test_that("no treatment_model with method='aiptw'", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
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
                                                         outcome_model=NULL,
                                                         cause=1),
               NULL)
})

test_that("no outcome_model with method='aiptw'", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
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
                                                         treatment_model=NULL,
                                                         cause=1),
               NULL)
})

test_that("no outcome_model with method='direct'", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
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
                                                         cause=1),
               NULL)
})

outcome_model <- list(cause=2)
class(outcome_model) <- "FGR"

test_that("cause is not the same as the cause in FGR", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
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
                                                    outcome_model=outcome_model,
                                                    cause=1),
               NULL)
})

test_that("conf_int with non-supported method", {
  expect_warning(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                         variable="group",
                                                         ev_time="time",
                                                         event="event",
                                                         method="matching",
                                                         conf_int=TRUE,
                                                         conf_level=0.95,
                                                         times=NULL,
                                                         bootstrap=FALSE,
                                                         n_boot=2,
                                                         na.action="na.omit",
                                                         clean_data=TRUE,
                                                         treatment_model=NULL,
                                                         cause=1),
               NULL)
})

## multiple imputation stuff

sim_dat_na_time <- sim_dat_na_group <- sim_dat_na_event <- sim_dat_na_x1 <-
  sim_dat

sim_dat_na_time$time[3] <- NA
sim_dat_na_group$group[3] <- NA
sim_dat_na_event$event[3] <- NA
sim_dat_na_x1$x1[3] <- NA

mids_na_time <- mice::mice(sim_dat_na_time, m=2, printFlag=FALSE)
mids_na_group <- mice::mice(sim_dat_na_group, m=2, printFlag=FALSE)
mids_na_event <- mice::mice(sim_dat_na_event, m=2, printFlag=FALSE)
mids_na_x1 <- mice::mice(sim_dat_na_x1, m=2, printFlag=FALSE)

outcome_model <- list()
class(outcome_model) <- "CauseSpecificCox"

censoring_model <- list()
class(censoring_model) <- "coxph"

treatment_model <- list()
class(treatment_model) <- "glm"

test_that("wrong outcome_model with mira", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=mids_na_x1,
                                                   variable="group",
                                                   ev_time="time",
                                                   event="event",
                                                   method="direct",
                                                   conf_int=TRUE,
                                                   conf_level=0.95,
                                                   times=NULL,
                                                   bootstrap=FALSE,
                                                   n_boot=2,
                                                   cause=1,
                                                   na.action="na.omit",
                                                   clean_data=TRUE,
                                                   outcome_model=outcome_model),
               NULL)
})

test_that("wrong treatment_model with mira", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=mids_na_x1,
                                              variable="group",
                                              ev_time="time",
                                              event="event",
                                              method="iptw",
                                              conf_int=TRUE,
                                              conf_level=0.95,
                                              times=NULL,
                                              bootstrap=FALSE,
                                              n_boot=2,
                                              na.action="na.omit",
                                              clean_data=TRUE,
                                              cause=1,
                                              treatment_model=treatment_model),
               NULL)
})

test_that("wrong censoring_model with mira", {
  expect_error(adjustedCurves:::check_inputs_adjustedcif(data=mids_na_x1,
                                               variable="group",
                                               ev_time="time",
                                               event="event",
                                               method="aiptw",
                                               conf_int=TRUE,
                                               conf_level=0.95,
                                               times=NULL,
                                               bootstrap=FALSE,
                                               n_boot=2,
                                               cause=1,
                                               na.action="na.omit",
                                               clean_data=TRUE,
                                               censoring_model=censoring_model),
               NULL)
})

test_that("warning with missing values in variable", {
  expect_warning(adjustedCurves:::check_inputs_adjustedcif(data=mids_na_group,
                                                        variable="group",
                                                        ev_time="time",
                                                        event="event",
                                                        method="aalen_johansen",
                                                        conf_int=TRUE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=FALSE,
                                                        n_boot=2,
                                                        cause=1,
                                                        na.action="na.omit",
                                                        clean_data=TRUE),
                 NULL)
})

test_that("warning with missing values in ev_time", {
  expect_warning(adjustedCurves:::check_inputs_adjustedcif(
    data=mids_na_time,
    variable="group",
    ev_time="time",
    event="event",
    method="aalen_johansen",
    conf_int=TRUE,
    conf_level=0.95,
    times=NULL,
    bootstrap=FALSE,
    n_boot=2,
    cause=1,
    na.action="na.omit",
    clean_data=TRUE),
    NULL)
})

test_that("warning with missing values in event", {
  expect_warning(adjustedCurves:::check_inputs_adjustedcif(data=mids_na_event,
                                                        variable="group",
                                                        ev_time="time",
                                                        event="event",
                                                        method="aalen_johansen",
                                                        conf_int=TRUE,
                                                        conf_level=0.95,
                                                        times=NULL,
                                                        bootstrap=FALSE,
                                                        n_boot=2,
                                                        cause=1,
                                                        clean_data=TRUE,
                                                        na.action="na.omit"),
                 NULL)
})
