
suppressMessages(requireNamespace("survival"))

set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_surv_n_100.Rds",
                               package="adjustedCurves"))
sim_dat$group <- factor(sim_dat$group)

adj <- adjustedsurv(data=sim_dat,
                    variable="group",
                    ev_time="time",
                    event="event",
                    method="km",
                    conf_int=TRUE,
                    bootstrap=TRUE,
                    n_boot=2)

adj2 <- adj
adj2$adjcif <- adj$adjsurv
class(adj2) <- "adjustedcif"

test_that("not an adjustedsurv/adjustedcif object", {
  expect_error(check_inputs_plot_difference(x="a",
                                            group_1="0",
                                            group_2="1",
                                            conf_int=TRUE,
                                            type="lines",
                                            max_t=Inf),
               paste0("'x' must be an 'adjustedsurv' object created using",
                      " the adjustedsurv function or an 'adjustedcif' object ",
                      "created using the adjustedcif function."),
               fixed=TRUE)
})

test_that("wrong group_1", {
  expect_error(check_inputs_plot_difference(x=adj,
                                            group_1=c("0", "1"),
                                            group_2="1",
                                            conf_int=TRUE,
                                            type="lines",
                                            max_t=Inf),
               paste0("'group_1' has to be a single character vector, ",
                      "specifying one of the treatment groups in 'variable'."),
               fixed=TRUE)
})

test_that("wrong group_2", {
  expect_error(check_inputs_plot_difference(x=adj,
                                            group_1="0",
                                            group_2=c("0", "1"),
                                            conf_int=TRUE,
                                            type="lines",
                                            max_t=Inf),
               paste0("'group_2' has to be a single character vector, ",
                      "specifying one of the treatment groups in 'variable'."),
               fixed=TRUE)
})

test_that("group_1 not in adjustedsurv", {
  expect_error(check_inputs_plot_difference(x=adj,
                                            group_1="2",
                                            group_2="0",
                                            conf_int=TRUE,
                                            type="lines",
                                            max_t=Inf),
               paste0("2 is not a valid group in 'variable'."),
               fixed=TRUE)
})

test_that("group_1 not in adjustedcif", {
  expect_error(check_inputs_plot_difference(x=adj2,
                                            group_1="2",
                                            group_2="0",
                                            conf_int=TRUE,
                                            type="lines",
                                            max_t=Inf),
               paste0("2 is not a valid group in 'variable'."),
               fixed=TRUE)
})

test_that("group_2 not in adjustedsurv", {
  expect_error(check_inputs_plot_difference(x=adj,
                                            group_1="0",
                                            group_2="2",
                                            conf_int=TRUE,
                                            type="lines",
                                            max_t=Inf),
               paste0("2 is not a valid group in 'variable'."),
               fixed=TRUE)
})

test_that("group_2 not in adjustedcif", {
  expect_error(check_inputs_plot_difference(x=adj2,
                                            group_1="0",
                                            group_2="2",
                                            conf_int=TRUE,
                                            type="lines",
                                            max_t=Inf),
               paste0("2 is not a valid group in 'variable'."),
               fixed=TRUE)
})

test_that("group_1 same as group_2", {
  expect_error(check_inputs_plot_difference(x=adj,
                                            group_1="0",
                                            group_2="0",
                                            conf_int=TRUE,
                                            type="lines",
                                            max_t=Inf),
               paste0("'group_1' and 'group_2' may not be equal."),
               fixed=TRUE)
})

test_that("undefined type", {
  expect_error(check_inputs_plot_difference(x=adj,
                                            group_1="0",
                                            group_2="1",
                                            conf_int=TRUE,
                                            type="undefined",
                                            max_t=Inf),
               paste0("'type' must be a single character string equal to one",
                      " of:  c('steps', 'lines', 'points', 'none')."),
               fixed=TRUE)
})

test_that("wrong max_t", {
  expect_error(check_inputs_plot_difference(x=adj,
                                            group_1="0",
                                            group_2="1",
                                            conf_int=TRUE,
                                            type="lines",
                                            max_t=-2),
               paste0("'max_t' must be a single number bigger than 0."),
               fixed=TRUE)
})

adj$boot_data <- NULL

test_that("no conf_int without bootstrapping", {
  expect_error(check_inputs_plot_difference(x=adj,
                                            group_1="0",
                                            group_2="1",
                                            conf_int=TRUE,
                                            type="lines",
                                            max_t=Inf),
               paste0("Confidence interval calculation is only possible ",
                      "when bootstrap=TRUE was used in the original ",
                      "adjustedsurv or adjustedcif function call."),
               fixed=TRUE)
})

adj$mids_analyses <- 1

test_that("no conf_int with MI", {
  expect_error(check_inputs_plot_difference(x=adj,
                                            group_1="0",
                                            group_2="1",
                                            conf_int=TRUE,
                                            type="lines",
                                            max_t=Inf),
               paste0("Confidence interval calculation is currently not ",
                      "possible when multiple imputation was used."),
               fixed=TRUE)
})
