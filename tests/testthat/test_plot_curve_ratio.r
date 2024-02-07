
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

test_that("plot, no arguments", {
  plt <- plot_curve_ratio(adj, group_1="1", group_2="0")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, no arguments", fig=plt)
})

test_that("plot, reversed differences", {
  plt <- plot_curve_ratio(adj)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, reversed differences", fig=plt)
})

test_that("plot, with conf_int", {
  plt <- suppressWarnings(plot_curve_ratio(adj, conf_int=TRUE))
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with conf_int", fig=plt)
})

test_that("plot, with conf_int boot", {
  plt <- plot_curve_ratio(adj, conf_int=TRUE, use_boot=TRUE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with conf_int boot", fig=plt)
})

test_that("plot, with lines", {
  plt <- plot_curve_ratio(adj, type="lines")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with lines", fig=plt)
})

test_that("plot, with lines ci", {
  plt <- suppressWarnings(plot_curve_ratio(adj, type="lines", conf_int=TRUE))
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with lines ci", fig=plt)
})

test_that("plot, without line at 1", {
  plt <- plot_curve_ratio(adj, line_at_ref=FALSE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, without line at ref", fig=plt)
})

sim_dat$event[1] <- 2
adj <- adjustedcif(data=sim_dat,
                   variable="group",
                   ev_time="time",
                   event="event",
                   method="aalen_johansen",
                   conf_int=TRUE,
                   bootstrap=TRUE,
                   n_boot=2,
                   cause=1)

test_that("plot, cif no arguments", {
  plt <- plot_curve_ratio(adj, group_1="1", group_2="0")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, cif no arguments", fig=plt)
})
