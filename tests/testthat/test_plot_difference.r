
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
  plt <- plot_difference(adj, group_1="1", group_2="0")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, no arguments", fig=plt)
})

test_that("plot, reversed differences", {
  plt <- plot_difference(adj)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, reversed differences", fig=plt)
})

test_that("plot, with conf_int", {
  plt <- plot_difference(adj, conf_int=TRUE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with conf_int", fig=plt)
})

test_that("plot, with lines", {
  plt <- plot_difference(adj, type="lines")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with lines", fig=plt)
})

test_that("plot, with lines ci", {
  plt <- plot_difference(adj, type="lines", conf_int=TRUE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with lines ci", fig=plt)
})

test_that("plot, with points", {
  plt <- plot_difference(adj, type="points")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with points", fig=plt)
})

test_that("plot, with points ci", {
  plt <- plot_difference(adj, type="points", conf_int=TRUE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with points ci", fig=plt)
})

test_that("plot, with none", {
  plt <- plot_difference(adj, type="none")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with none", fig=plt)
})

test_that("plot, with loess", {
  plt <- plot_difference(adj, type="none", loess_smoother=TRUE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with loess", fig=plt)
})

test_that("plot, without line at 0", {
  plt <- plot_difference(adj, line_at_0=FALSE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, without line at 0", fig=plt)
})

test_that("plot, with fill_area lines", {
  plt <- plot_difference(adj, type="lines", fill_area=TRUE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with fill_area lines", fig=plt)
})

test_that("plot, with fill_area steps", {
  plt <- plot_difference(adj, type="steps", fill_area=TRUE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with fill_area steps", fig=plt)
})

test_that("plot, with much stuff", {
  plt <- plot_difference(adj, conf_int=TRUE, color="blue", linetype="dotted",
                         alpha=0.8, line_at_0_size=1.1, line_at_0_color="red",
                         loess_smoother=TRUE, loess_span=0.55)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with much stuff", fig=plt)
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
  plt <- plot_difference(adj, group_1="1", group_2="0")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, cif no arguments", fig=plt)
})
