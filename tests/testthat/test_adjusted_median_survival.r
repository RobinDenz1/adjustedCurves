
suppressMessages(requireNamespace("survival"))

set.seed(42)
sim_dat <- readRDS(system.file("testdata",
                               "d_sim_surv_n_50.Rds",
                               package="adjustedCurves"))
sim_dat$group <- factor(sim_dat$group)

adj <- adjustedsurv(data=sim_dat,
                    variable="group",
                    ev_time="time",
                    event="event",
                    method="km",
                    conf_int=TRUE,
                    bootstrap=FALSE)

test_that("median surv, no boot", {
  adj_med <- adjusted_median_survival(adj, verbose=FALSE)
  expect_equal(round(adj_med$median_surv, 4), c(0.4785, 0.6252))
})

test_that("median surv, with verbose", {
  expect_snapshot_output(adjusted_median_survival(adj, verbose=TRUE))
})
