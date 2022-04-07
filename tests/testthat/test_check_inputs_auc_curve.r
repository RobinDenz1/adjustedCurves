
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

test_that("wrong times", {
  expect_error(check_inputs_auc_curve(times=list(1, 2),
                                      max_t=Inf,
                                      color=TRUE,
                                      linetype=TRUE,
                                      facet=TRUE),
               paste0("'times' must be a numeric vector containing only ",
                      "values > 0 or NULL."),
               fixed=TRUE)
})

test_that("wrong max_t", {
  expect_error(check_inputs_auc_curve(times=NULL,
                                      max_t=0,
                                      color=TRUE,
                                      linetype=TRUE,
                                      facet=TRUE),
               "'max_t' must be a single number > 0.",
               fixed=TRUE)
})

test_that("wrong color", {
  expect_error(check_inputs_auc_curve(times=NULL,
                                      max_t=Inf,
                                      color="black",
                                      linetype=TRUE,
                                      facet=TRUE),
               paste0("'color' must be either TRUE or FALSE. To use custom",
                      " colors, the 'custom_colors' argument should be used."),
               fixed=TRUE)
})

test_that("wrong linetype", {
  expect_error(check_inputs_auc_curve(times=NULL,
                                      max_t=Inf,
                                      color=TRUE,
                                      linetype="dashed",
                                      facet=TRUE),
               paste0("'linetype' must be either TRUE or FALSE. To use custom",
                      " linetypes, the 'custom_linetypes' argument should",
                      " be used."),
               fixed=TRUE)
})

test_that("wrong facet", {
  expect_error(check_inputs_auc_curve(times=NULL,
                                      max_t=Inf,
                                      color=FALSE,
                                      linetype=TRUE,
                                      facet="A"),
               "'facet' must be either TRUE or FALSE.",
               fixed=TRUE)
})
