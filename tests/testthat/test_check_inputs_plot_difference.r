
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

adj_test_cat <- list(categorical=TRUE)
class(adj_test_cat) <- "curve_test"

test_that("not an adjustedsurv/adjustedcif object", {
  expect_error(check_inputs_plot_difference(x="a",
                                            group_1="0",
                                            group_2="1",
                                            conf_int=TRUE,
                                            type="lines",
                                            max_t=Inf,
                                            use_boot=FALSE,
                                            test=NULL,
                                            integral_from=0,
                                            integral_to=NULL,
                                            p_value=FALSE,
                                            integral=FALSE),
               paste0("'x' must be an 'adjustedsurv' object created using",
                      " the adjustedsurv function or an 'adjustedcif' object ",
                      "created using the adjustedcif function."),
               fixed=TRUE)
})

test_that("undefined type", {
  expect_error(check_inputs_plot_difference(x=adj,
                                            group_1="0",
                                            group_2="1",
                                            conf_int=TRUE,
                                            type="undefined",
                                            max_t=Inf,
                                            use_boot=FALSE,
                                            test=NULL,
                                            integral_from=0,
                                            integral_to=NULL,
                                            p_value=FALSE,
                                            integral=FALSE),
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
                                            max_t=-2,
                                            use_boot=FALSE,
                                            test=NULL,
                                            integral_from=0,
                                            integral_to=NULL,
                                            p_value=FALSE,
                                            integral=FALSE),
               paste0("'max_t' must be a single number bigger than 0."),
               fixed=TRUE)
})

test_that("wrong test", {
  expect_error(check_inputs_plot_difference(x=adj,
                                            group_1="0",
                                            group_2="1",
                                            conf_int=TRUE,
                                            type="lines",
                                            max_t=Inf,
                                            use_boot=FALSE,
                                            test="",
                                            integral_from=0,
                                            integral_to=NULL,
                                            p_value=FALSE,
                                            integral=FALSE),
               paste0("'test' must be either NULL or a 'curve_test' object ",
                      "created using the adjusted_curve_test function."),
               fixed=TRUE)
})

test_that("test with categorical", {
  expect_error(check_inputs_plot_difference(x=adj,
                                            group_1="0",
                                            group_2="1",
                                            conf_int=TRUE,
                                            type="lines",
                                            max_t=Inf,
                                            use_boot=FALSE,
                                            test=adj_test_cat,
                                            integral_from=0,
                                            integral_to=NULL,
                                            p_value=FALSE,
                                            integral=FALSE),
               paste0("The curve_test object supplied in the 'test' ",
                      "argument may only contain two groups, corresponding ",
                      "to the groups used in the plot_difference function."),
               fixed=TRUE)
})

test_that("wrong integral_from", {
  expect_error(check_inputs_plot_difference(x=adj,
                                            group_1="0",
                                            group_2="1",
                                            conf_int=TRUE,
                                            type="lines",
                                            max_t=Inf,
                                            use_boot=FALSE,
                                            test=NULL,
                                            integral_from=-1,
                                            integral_to=NULL,
                                            p_value=FALSE,
                                            integral=FALSE),
               paste0("'integral_from' must be a single number > 0 or NULL."),
               fixed=TRUE)
})

test_that("wrong integral_to", {
  expect_error(check_inputs_plot_difference(x=adj,
                                            group_1="0",
                                            group_2="1",
                                            conf_int=TRUE,
                                            type="lines",
                                            max_t=Inf,
                                            use_boot=FALSE,
                                            test=NULL,
                                            integral_from=0,
                                            integral_to="A",
                                            p_value=FALSE,
                                            integral=FALSE),
               paste0("'integral_to' must be a single number."),
               fixed=TRUE)
})

test_that("integral_to not smaller than integral_from", {
  expect_error(check_inputs_plot_difference(x=adj,
                                            group_1="0",
                                            group_2="1",
                                            conf_int=TRUE,
                                            type="lines",
                                            max_t=Inf,
                                            use_boot=FALSE,
                                            test=NULL,
                                            integral_from=2.1,
                                            integral_to=2,
                                            p_value=FALSE,
                                            integral=FALSE),
               paste0("'integral_from' must be smaller than 'integral_to'."),
               fixed=TRUE)
})

test_that("integral without needed args", {
  expect_error(check_inputs_plot_difference(x=adj,
                                            group_1="0",
                                            group_2="1",
                                            conf_int=TRUE,
                                            type="lines",
                                            max_t=Inf,
                                            use_boot=FALSE,
                                            test=NULL,
                                            integral_from=0,
                                            integral_to=NULL,
                                            p_value=FALSE,
                                            integral=TRUE),
               paste0("If 'integral' is specified, either 'test' or ",
                      "'integral_to' also need to be specified. See details."),
               fixed=TRUE)
})

test_that("p_value without needed args", {
  expect_error(check_inputs_plot_difference(x=adj,
                                            group_1="0",
                                            group_2="1",
                                            conf_int=TRUE,
                                            type="lines",
                                            max_t=Inf,
                                            use_boot=FALSE,
                                            test=NULL,
                                            integral_from=0,
                                            integral_to=NULL,
                                            p_value=TRUE,
                                            integral=FALSE),
               paste0("If 'p_value' is specified, either 'test' or ",
                      "'integral_to' also need to be specified. See details."),
               fixed=TRUE)
})

adj$boot_adjsurv <- NULL

test_that("using p_value without boot", {
  expect_error(check_inputs_plot_difference(x=adj,
                                            group_1="0",
                                            group_2="1",
                                            conf_int=FALSE,
                                            type="lines",
                                            max_t=Inf,
                                            use_boot=FALSE,
                                            test=NULL,
                                            integral_from=0,
                                            integral_to=1,
                                            p_value=TRUE,
                                            integral=FALSE),
               paste0("'p_value' can only be used when bootstrap=TRUE was ",
                      "used in the original adjustedsurv or adjustedcif ",
                      "function call."),
               fixed=TRUE)
})

test_that("using use_boot without boot", {
  expect_error(check_inputs_plot_difference(x=adj,
                                            group_1="0",
                                            group_2="1",
                                            conf_int=TRUE,
                                            type="lines",
                                            max_t=Inf,
                                            use_boot=TRUE,
                                            test=NULL,
                                            integral_from=0,
                                            integral_to=NULL,
                                            p_value=FALSE,
                                            integral=FALSE),
               paste0("Bootstrapped estimates can only be calculated if ",
                      "'bootstrap=TRUE' was used in the original ",
                      "adjustedsurv or adjustedcif function call."),
               fixed=TRUE)
})

adj$adjsurv <- dplyr::select(adj$adjsurv, c("time", "group", "surv"))

test_that("using conf_int without approximate stuff surv", {
  expect_error(check_inputs_plot_difference(x=adj,
                                            group_1="0",
                                            group_2="1",
                                            conf_int=TRUE,
                                            type="lines",
                                            max_t=Inf,
                                            use_boot=FALSE,
                                            test=NULL,
                                            integral_from=0,
                                            integral_to=NULL,
                                            p_value=FALSE,
                                            integral=FALSE),
               paste0("There are no approximate standard error ",
                      "calculations to use. Either set 'use_boot=TRUE' or ",
                      "rerun the adjustedsurv function with ",
                      "'conf_int=TRUE' if possible."),
               fixed=TRUE)
})

adj2$adjcif <- dplyr::select(adj2$adjcif, c("time", "group", "surv"))

test_that("using conf_int without approximate stuff surv", {
  expect_error(check_inputs_plot_difference(x=adj2,
                                            group_1="0",
                                            group_2="1",
                                            conf_int=TRUE,
                                            type="lines",
                                            max_t=Inf,
                                            use_boot=FALSE,
                                            test=NULL,
                                            integral_from=0,
                                            integral_to=NULL,
                                            p_value=FALSE,
                                            integral=FALSE),
               paste0("There are no approximate standard error ",
                      "calculations to use. Either set 'use_boot=TRUE' or ",
                      "rerun the adjustedcif function with ",
                      "'conf_int=TRUE' if possible."),
               fixed=TRUE)
})
