
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
  plt <- plot_curve_diff(adj, group_1="1", group_2="0")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, no arguments", fig=plt)
})

test_that("plot, reversed differences", {
  plt <- plot_curve_diff(adj)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, reversed differences", fig=plt)
})

test_that("plot, with conf_int", {
  plt <- plot_curve_diff(adj, conf_int=TRUE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with conf_int", fig=plt)
})

test_that("plot, with conf_int boot", {
  plt <- plot_curve_diff(adj, conf_int=TRUE, use_boot=TRUE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with conf_int boot", fig=plt)
})

test_that("plot, with lines", {
  plt <- plot_curve_diff(adj, type="lines")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with lines", fig=plt)
})

test_that("plot, with lines ci", {
  plt <- plot_curve_diff(adj, type="lines", conf_int=TRUE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with lines ci", fig=plt)
})

test_that("plot, with points", {
  plt <- plot_curve_diff(adj, type="points")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with points", fig=plt)
})

test_that("plot, with points ci", {
  plt <- plot_curve_diff(adj, type="points", conf_int=TRUE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with points ci", fig=plt)
})

test_that("plot, with none", {
  plt <- plot_curve_diff(adj, type="none")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with none", fig=plt)
})

test_that("plot, with loess", {
  plt <- plot_curve_diff(adj, type="none", loess_smoother=TRUE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with loess", fig=plt)
})

test_that("plot, without line at 0", {
  plt <- plot_curve_diff(adj, line_at_0=FALSE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, without line at 0", fig=plt)
})

test_that("plot, with fill_area lines", {
  plt <- plot_curve_diff(adj, type="lines", fill_area=TRUE,
                         fill_only_interval=FALSE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with fill_area lines", fig=plt)
})

test_that("plot, with fill_area steps", {
  plt <- plot_curve_diff(adj, type="steps", fill_area=TRUE,
                         fill_only_interval=FALSE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with fill_area steps", fig=plt)
})

test_that("plot, with p_value", {
  plt <- plot_curve_diff(adj, type="steps", p_value=TRUE, integral_to=0.3)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with p_value", fig=plt)
})

test_that("plot, with p_value linear", {
  plt <- plot_curve_diff(adj, type="lines", p_value=TRUE, integral_to=0.3)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with p_value linear", fig=plt)
})

test_that("plot, with integral", {
  plt <- plot_curve_diff(adj, type="steps", integral=TRUE, integral_to=0.3)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with integral", fig=plt)
})

test_that("plot, with integral linear", {
  plt <- plot_curve_diff(adj, type="lines", integral=TRUE, integral_to=0.3)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with integral linear", fig=plt)
})

test_that("plot, with all texts", {
  plt <- plot_curve_diff(adj, type="steps", p_value=TRUE, integral=TRUE,
                         integral_to=0.3, interval=TRUE, text_pos_x="right",
                         text_pos_y="top", text_size=4, text_family="serif",
                         text_fontface="bold", text_color="blue",
                         text_alpha=0.6, text_digits=5, text_format_p=FALSE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with all texts", fig=plt)
})

adj_test <- adjusted_curve_test(adj, to=0.75, from=0, interpolation="steps")

test_that("plot, with p_value test", {
  plt <- plot_curve_diff(adj, type="steps", p_value=TRUE, test=adj_test)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with p_value test", fig=plt)
})

test_that("plot, with integral test", {
  plt <- plot_curve_diff(adj, type="steps", integral=TRUE, test=adj_test)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with integral test", fig=plt)
})

test_that("plot, with much stuff", {
  plt <- plot_curve_diff(adj, conf_int=TRUE, color="blue", linetype="dotted",
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
  plt <- plot_curve_diff(adj, group_1="1", group_2="0")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, cif no arguments", fig=plt)
})
