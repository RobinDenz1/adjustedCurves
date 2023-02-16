
suppressMessages(requireNamespace("survival"))

set.seed(42)
sim_dat <- readRDS(system.file("testdata",
                               "d_sim_crisk_n_20.Rds",
                               package="adjustedCurves"))
sim_dat$group <- as.factor(sim_dat$group)

test_that("variable not in data", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               "grouup is not a valid column name in 'data'.")
})


test_that("variable has wrong type", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               paste0("The column in 'data' specified by 'variable' ",
                      "needs to be a factor variable."), fixed=TRUE)
})

test_that("variable has really wrong type", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
                                        variable=1,
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
               paste0("Arguments 'variable', 'ev_time', 'event' and ",
                      "'method' must be character strings, specifying ",
                      "variables in 'data'."), fixed=TRUE)
})

test_that("variable has wrong type with matching", {
  sim_dat_err <- sim_dat
  sim_dat_err$group <- as.character(sim_dat_err$group)
  expect_error(check_inputs_adjustedcif(data=sim_dat_err,
                                        variable="group",
                                        ev_time="time",
                                        event="event",
                                        method="matching",
                                        conf_int=FALSE,
                                        conf_level=0.95,
                                        times=NULL,
                                        bootstrap=FALSE,
                                        n_boot=2,
                                        na.action="na.omit",
                                        clean_data=TRUE,
                                        cause=1),
               paste0("The column in 'data' specified by 'variable' needs ",
                      "to be a factor variable."),
               fixed=TRUE)
})

test_that("variable has wrong length", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               paste0("'variable' must be a character string of length 1 ",
                      "specifying the grouping variable in 'data'."),
               fixed=TRUE)
})

test_that("non-standard evaluation in variable", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               paste0("'variable' must be a character string ",
                      "specifying a variable in 'data'."))
})

test_that("variable is a constant", {
  sim_dat_err <- sim_dat
  sim_dat_err$group <- factor(1)
  expect_error(check_inputs_adjustedcif(data=sim_dat_err,
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
                                        cause=1),
               "There have to be at least two groups in 'variable'.",
               fixed=TRUE)
})

test_that("categorical variable not allowed", {
  sim_dat_err <- sim_dat
  sim_dat_err$group <- c(rep(0, 10), rep(1, 5), rep(2, 5))
  sim_dat_err$group <- as.factor(sim_dat_err$group)
  expect_error(check_inputs_adjustedcif(data=sim_dat_err,
                                        variable="group",
                                        ev_time="time",
                                        event="event",
                                        method="matching",
                                        conf_int=FALSE,
                                        conf_level=0.95,
                                        times=NULL,
                                        bootstrap=FALSE,
                                        n_boot=2,
                                        na.action="na.omit",
                                        clean_data=TRUE,
                                        cause=1),
  "Categorical treatments are currently not supported for method='matching'.")
})

test_that("ev_time not in data", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               paste0("time2 is not a valid column name in 'data'."))
})

test_that("ev_time wrong format", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               "The column in 'data' specified by 'ev_time' must be numeric.")
})

sim_dat$time2 <- sim_dat$time
sim_dat$time2[1] <- -1

test_that("ev_time includes negative values", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               paste0("'ev_time' must be a character string of length 1 ",
                      "specifying the time until an event or ",
                      "censoring occured."))
})

test_that("non-standard evaluation in ev_time", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
                                        variable="group",
                                        ev_time=event_time,
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
               paste0("'ev_time' must be a character string specifying ",
                      "a variable in 'data'."))
})

test_that("event not in data", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               "event2 is not a valid column name in 'data'.")
})

test_that("event wrong format", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               "The column in 'data' specified by 'event' must be numeric.")
})

sim_dat$event2 <- sample(c(0, 1), size=nrow(sim_dat), replace=TRUE)

test_that("event only includes two values", {
  expect_warning(check_inputs_adjustedcif(data=sim_dat,
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
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               paste0("'event' must be a character string of length 1 ",
                      "specifying the numeric event indicator."))
})

test_that("non-standard evaluation in event", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
                                        variable="group",
                                        ev_time="time",
                                        event=status,
                                        method="aalen_johansen",
                                        conf_int=FALSE,
                                        conf_level=0.95,
                                        times=NULL,
                                        bootstrap=FALSE,
                                        n_boot=2,
                                        na.action="na.omit",
                                        clean_data=TRUE,
                                        cause=1),
               paste0("'event' must be a character string ",
                      "specifying a variable in 'data'."))
})

test_that("cause has wrong type", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               paste0("'cause' must be a number specifying the cause of ",
                      "interest in the column specified with 'event'."))
})

test_that("more than one cause supplied", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               "'cause' must be of length = 1.")
})

test_that("method undefined", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               paste0("Method 'km3' is undefined. See documentation ",
                      "for details on available methods."))
})

test_that("wrong method length", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               paste0("'method' must be a single character string. Using ",
                      "multiple methods in one call is currently ",
                      "not supported."))
})

test_that("wrong conf_int type", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               paste0("'conf_int' must be either TRUE or FALSE."))
})

test_that("wrong conf_int length", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
              paste0("'conf_int' must be either TRUE or FALSE, not a vector."))
})

test_that("wrong conf_level type", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
                                        variable="group",
                                        ev_time="time",
                                        event="event",
                                        method="aalen_johansen",
                                        conf_int=TRUE,
                                        conf_level="0.95",
                                        times=NULL,
                                        bootstrap=FALSE,
                                        n_boot=2,
                                        na.action="na.omit",
                                        clean_data=TRUE,
                                        cause=1),
               "'conf_level' must be a number < 1 and > 0.")
})

test_that("conf_level not valid", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               "'conf_level' must be a number < 1 and > 0.")
})

test_that("wrong conf_level length", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               "'conf_level' must be a number < 1 and > 0, not a vector.")
})

test_that("wrong bootstrap type", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               "'bootstrap' must be either TRUE or FALSE.")
})

test_that("wrong bootstrap length", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               "'bootstrap' must be either TRUE or FALSE, not a vector.")
})

test_that("wrong n_boot type", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               "'n_boot' must be a positive integer > 1.")
})

test_that("wrong n_boot length", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               "'n_boot' must be a positive integer > 2, not a vector.")
})

test_that("wrong times", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               "'times' must be a numeric vector or NULL.")
})

test_that("wrong na.action type", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
    "'na.action' must be a function or a single character string. See details.")
})

test_that("wrong na.action length", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               paste0("'na.action' must be a function or a character ",
                      "string, not a vector. See documentation."))
})

test_that("wrong clean_data type", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               "'clean_data' must be either TRUE or FALSE.")
})

test_that("wrong clean_data length", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               "'clean_data' must be either TRUE or FALSE, not a vector.")
})

test_that("no extrapolation", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               paste0("Values in 'times' must be smaller than ",
                      "max(data[,ev_time]). No extrapolation allowed."),
               fixed=TRUE)
})

test_that("bootstrap only with treatment_model that can be updated", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               paste0("'treatment_model' needs to be a model that can ",
                      "be refit or a formula object when using ",
                      "bootstrap=TRUE."))
})

test_that("no outcome_vars argument", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
                                        cause=1),
               paste0("Argument 'outcome_vars' must be specified when using ",
                      "method='direct_pseudo'. See documentation."))
})

test_that("wrong outcome_vars argument", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               paste0("'outcome_vars' should be a character vector of ",
                      "column names in 'data', used to model ",
                      "the outcome mechanism."))
})

test_that("wrong type_time argument", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               "'type_time' should be either 'factor', 'bs' or 'ns'.")
})

test_that("too many spline_df", {
  expect_warning(check_inputs_adjustedcif(data=sim_dat,
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
                 "'spline_df' > len(times) might lead to problems.",
                 fixed=TRUE)
})

test_that("too many spline_df default", {
  expect_warning(check_inputs_adjustedcif(data=sim_dat,
                                          variable="group",
                                          ev_time="time",
                                          event="event",
                                          method="direct_pseudo",
                                          conf_int=FALSE,
                                          conf_level=0.95,
                                          bootstrap=TRUE,
                                          n_boot=2,
                                          na.action="na.omit",
                                          clean_data=TRUE,
                                          outcome_vars=c("x1"),
                                          type_time="bs",
                                          times=c(0.1, 0.2, 0.4),
                                          cause=1),
                 paste0("'spline_df' > length(times) might lead to problems ",
                        "when type_time!='factor'."),
                 fixed=TRUE)
})

test_that("not enough time points", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               paste0("'geese' models require at least two distinct ",
                      "time points. Add more points in time to ",
                      "'times' and run again."))
})

test_that("weights in matching", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               paste0("Propensity Scores > 1 or < 0 not allowed. ",
                      "Perhaps you supplied weights on accident?"))
})

test_that("ps_scores > 1 or < 0", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
                                        variable="group",
                                        ev_time="time",
                                        event="event",
                                        method="aiptw_pseudo",
                                        conf_int=FALSE,
                                        conf_level=0.95,
                                        times=NULL,
                                        bootstrap=TRUE,
                                        n_boot=2,
                                        na.action="na.omit",
                                        clean_data=TRUE,
                                        outcome_vars=c("x1"),
                                        cause=1,
                                        treatment_model=rep(2, nrow(sim_dat))),
               paste0("Propensity Scores supplied using the 'treatment_model' ",
                      "argument must be smaller than 1 and bigger than 0."))
})

test_that("ps_scores wrong length", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
                                        variable="group",
                                        ev_time="time",
                                        event="event",
                                        method="aiptw_pseudo",
                                        conf_int=FALSE,
                                        conf_level=0.95,
                                        times=NULL,
                                        bootstrap=TRUE,
                                        n_boot=2,
                                        na.action="na.omit",
                                        clean_data=TRUE,
                                        outcome_vars=c("x1"),
                                        cause=1,
                                        treatment_model=c(0.1, 0.2)),
               paste0("The vector of propensity score supplied in the ",
                      "'treatment_model' argument must be of ",
                      "length nrow(data)."), fixed=TRUE)
})

test_that("no treatment_model in matching", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               paste0("Argument 'treatment_model' must be specified ",
                      "when using method='matching'."))
})

test_that("no treatment_model with method='aiptw'", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
               paste0("Argument 'treatment_model' must be specified ",
                      "when using method='aiptw'."))
})

test_that("no 'treatment_model' with iptw", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
                                        variable="group",
                                        ev_time="time",
                                        event="event",
                                        method="iptw",
                                        conf_int=FALSE,
                                        conf_level=0.95,
                                        times=NULL,
                                        bootstrap=TRUE,
                                        n_boot=2,
                                        na.action="na.omit",
                                        clean_data=TRUE,
                                        cause=1),
               paste0("Argument 'treatment_model' must be defined when ",
                      "using method='iptw'. See documentation."))
})

test_that("wrong 'treatment_model' with iptw", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
                                        variable="group",
                                        ev_time="time",
                                        event="event",
                                        method="iptw",
                                        conf_int=FALSE,
                                        conf_level=0.95,
                                        times=NULL,
                                        bootstrap=TRUE,
                                        n_boot=2,
                                        na.action="na.omit",
                                        clean_data=TRUE,
                                        cause=1,
                                        treatment_model="haha"),
               paste0("'treatment_model' must be a glm or multinom object. ",
                      "See documentation."))
})

test_that("no outcome_model with method='aiptw'", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
      "Argument 'outcome_model' must be specified when using method='aiptw'.")
})

test_that("no 'treatment_model' with iptw_pseudo", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
                                        cause=1),
               paste0("Argument 'treatment_model' must be specified when ",
                      "using method='iptw_pseudo'. See documentation."))
})

test_that("no outcome_model with method='direct'", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
      "Argument 'outcome_model' must be specified when using method='direct'.")
})

test_that("wrong treatment_model with aiptw_pseudo", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
                                        variable="group",
                                        ev_time="time",
                                        event="event",
                                        method="aiptw_pseudo",
                                        conf_int=FALSE,
                                        conf_level=0.95,
                                        times=NULL,
                                        bootstrap=TRUE,
                                        n_boot=2,
                                        na.action="na.omit",
                                        clean_data=TRUE,
                                        cause=1,
                                        outcome_vars=c("x1", "x2", "x4"),
                                        treatment_model="wrong_treat_mod"),
               paste0("Argument 'treatment_model' must be one of: 'glm', ",
                      "'multinom', 'mira' or a numeric vector of ",
                      "propensity scores."))
})

outcome_model <- list(cause=2)
class(outcome_model) <- "FGR"

test_that("cause is not the same as the cause in FGR", {
  expect_error(check_inputs_adjustedcif(data=sim_dat,
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
              paste0("The FGR model needs to be fit with the same ",
                     "'cause' as specified in the 'cause' argument."))
})

test_that("conf_int with non-supported method", {
  expect_warning(check_inputs_adjustedcif(data=sim_dat,
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
               paste0("Asymptotic or exact variance calculations are ",
                      "currently not available for method='matching'. ",
                      "Use bootstrap=TRUE to get bootstrap estimates."))
})

## multiple imputation stuff

sim_dat_na_time <- sim_dat_na_group <- sim_dat_na_event <- sim_dat_na_x1 <-
  sim_dat

sim_dat_na_time$time[3] <- NA
sim_dat_na_group$group[3] <- NA
sim_dat_na_event$event[3] <- NA
sim_dat_na_x1$x1[3] <- NA

mids_na_time <- suppressWarnings(mice::mice(sim_dat_na_time, m=2,
                                            printFlag=FALSE))
mids_na_group <- suppressWarnings(mice::mice(sim_dat_na_group, m=2,
                                             printFlag=FALSE))
mids_na_event <- suppressWarnings(mice::mice(sim_dat_na_event, m=2,
                                             printFlag=FALSE))
mids_na_x1 <- suppressWarnings(mice::mice(sim_dat_na_x1, m=2, printFlag=FALSE))

outcome_model <- list()
class(outcome_model) <- "CauseSpecificCox"

censoring_model <- list()
class(censoring_model) <- "coxph"

treatment_model <- list()
class(treatment_model) <- "glm"

test_that("wrong outcome_model with mira", {
  expect_error(check_inputs_adjustedcif(data=mids_na_x1,
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
               paste0("When using multiple imputation, mira objects need ",
                      "to be supplied to 'outcome_model' ",
                      "instead of single models. See documentation."))
})

test_that("wrong treatment_model with mira", {
  expect_error(check_inputs_adjustedcif(data=mids_na_x1,
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
               paste0("When using multiple imputation, ",
                      "mira objects or a formula ",
                      "need to be supplied to 'treatment_model' instead of ",
                      "single models. See documentation."))
})

test_that("wrong censoring_model with mira", {
  expect_error(check_inputs_adjustedcif(data=mids_na_x1,
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
               paste0("When using multiple imputation, mira objects need to ",
                      "be supplied to 'censoring_model' instead of ",
                      "single models. See documentation."))
})

test_that("warning with missing values in variable", {
  expect_warning(check_inputs_adjustedcif(data=mids_na_group,
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
                 paste0("Using multiple imputation with missing values in ",
                        "'variable' has not been tested yet. ",
                        "Use with caution."))
})

test_that("warning with missing values in ev_time", {
  expect_warning(check_inputs_adjustedcif(data=mids_na_time,
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
    paste0("Using multiple imputation with missing values in 'ev_time' ",
           "variable has not been tested yet. Use with caution."))
})

test_that("warning with missing values in event", {
  expect_warning(check_inputs_adjustedcif(data=mids_na_event,
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
                 paste0("Using multiple imputation with missing values ",
                        "in 'event' variable has not been tested yet. ",
                        "Use with caution."))
})
