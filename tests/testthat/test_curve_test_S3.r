
set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_surv_n_20.Rds",
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

adj_test <- adjusted_curve_test(adj, to=0.2)

test_that("print.curve_test, default", {
  expect_snapshot_output(print(adj_test))
})

test_that("summary.curve_test, default", {
  expect_snapshot_output(summary(adj_test))
})

test_that("plot with wrong type", {
  adj_test$categorical <- TRUE
  expect_error(plot(adj_test, type="undefined"),
               "type='undefined' is not defined.",
               fixed=TRUE)
  adj_test$categorical <- FALSE
  expect_error(plot(adj_test, type="undefined"),
               "type='undefined' is not defined.",
               fixed=TRUE)
})

test_that("plot with mids not possible", {
  adj_test$mids_analyses <- list(A="I am not empty.")
  expect_error(plot(adj_test),
               paste0("There is no plot method for 'curve_test' objects ",
                      "fitted using multiply imputed datasets."),
               fixed=TRUE)
})
