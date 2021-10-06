
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
                                                         cause=1),
               NULL)
})

test_that("wrong conf_int", {
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
                                                        cause=1),
               NULL)
})

test_that("wrong conf_level", {
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
                                                        cause=1),
               NULL)
})

test_that("wrong bootstrap", {
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
                                                        cause=1),
               NULL)
})

test_that("wrong n_boot", {
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
                                                         outcome_vars=c("x1"),
                                                         cause=1),
               NULL)
})

test_that("bootstrap in matching", {
  expect_warning(adjustedCurves:::check_inputs_adjustedcif(data=sim_dat,
                                                           variable="group",
                                                           ev_time="time",
                                                           event="event",
                                                           method="matching",
                                                           conf_int=FALSE,
                                                           conf_level=0.95,
                                                           times=0.2,
                                                           bootstrap=TRUE,
                                                           n_boot=2,
                                                           na.action="na.omit",
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
                                      treatment_model=runif(20, min=10, max=20),
                                      cause=1),
               NULL)
})

test_that("no models with method='aiptw'", {
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
                                                         cause=1),
               NULL)
})

