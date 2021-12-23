library(survival)
library(vdiffr)

set.seed(42)

sim_dat <- sim_confounded_crisk(n=100)
sim_dat$group <- factor(sim_dat$group)

adj <- adjustedcif(data=sim_dat,
                   variable="group",
                   ev_time="time",
                   event="event",
                   method="aalen_johansen",
                   conf_int=TRUE,
                   bootstrap=TRUE,
                   n_boot=2,
                   cause=1)

test_that("plot, no arguments", {
  plt <- plot(adj)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, no arguments", fig=plt)
})

test_that("plot, with conf_int", {
  plt <- plot(adj, conf_int=TRUE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, with conf_int", fig=plt)
})

test_that("plot, using boot", {
  plt <- plot(adj, conf_int=TRUE, use_boot=TRUE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, using boot", fig=plt)
})

test_that("plot, using iso_reg", {
  plt <- plot(adj, iso_reg=TRUE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, using iso_reg", fig=plt)
})

test_that("plot, using censoring indicators (lines)", {
  plt <- plot(adj, iso_reg=TRUE,
              censoring_ind="lines",
              censoring_ind_width=0.1,
              censoring_ind_size=1,
              censoring_ind_alpha=0.5)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, using censoring indicators (lines)",
                              fig=plt)
})

test_that("plot, using censoring indicators (points)", {
  plt <- plot(adj, iso_reg=TRUE,
              censoring_ind="points",
              censoring_ind_size=1,
              censoring_ind_alpha=0.5,
              censoring_ind_shape=15)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, using censoring indicators (points)",
                              fig=plt)
})

test_that("plot, using max_t", {
  plt <- plot(adj, max_t=0.2)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, using max_t", fig=plt)
})

test_that("plot, using force_bounds", {
  plt <- plot(adj, force_bounds=TRUE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, using force_bounds", fig=plt)
})

test_that("plot, using linetype + facets + color", {
  plt <- plot(adj, color=TRUE, linetype=TRUE, facet=TRUE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, using linetype + facets + color", fig=plt)
})

test_that("plot, using labs + title", {
  plt <- plot(adj, xlab="X", ylab="Y", title="Title")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, using labs + title", fig=plt)
})

test_that("plot, using legend.title + legend.position", {
  plt <- plot(adj, legend.title="A", legend.position="bottom")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, using legend.title + legend.position",
                              fig=plt)
})

test_that("plot, using ggplot theme", {
  plt <- plot(adj, gg_theme=ggplot2::theme_bw())
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, using ggplot theme", fig=plt)
})

test_that("plot, using ylim", {
  plt <- plot(adj, ylim=c(0, 1))
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, using ylim", fig=plt)
})

test_that("plot, using single_color", {
  plt <- plot(adj, color=FALSE, single_color="blue")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, using single_color", fig=plt)
})

test_that("plot, using single_linetype", {
  plt <- plot(adj, linetype=FALSE, color=FALSE, single_linetype="dashed")
  expect_s3_class(plt, "ggplot")
})

test_that("plot, using custom_colors", {
  plt <- plot(adj, custom_colors=c("red", "blue"))
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, using custom_colors", fig=plt)
})

test_that("plot, using custom_linetypes", {
  plt <- plot(adj, custom_linetypes=c("dashed", "solid"), linetype=TRUE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, using custom_linetypes", fig=plt)
})

test_that("plot, using ci_draw_alpha", {
  plt <- plot(adj, conf_int=TRUE, ci_draw_alpha=0.1)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, using ci_draw_alpha", fig=plt)
})

test_that("plot, using steps", {
  plt <- plot(adj, steps=FALSE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, using steps", fig=plt)
})

test_that("plot, using many many things", {
  plt <- plot(adj,
              conf_int=TRUE,
              max_t=0.8,
              use_boot=FALSE,
              force_bounds=TRUE,
              iso_reg=TRUE,
              color=TRUE,
              linetype=TRUE,
              facet=TRUE,
              line_size=1.2,
              xlab="X",
              ylab="Y",
              title="Title",
              legend.title="Legend Title",
              legend.position="bottom",
              gg_theme=ggplot2::theme_bw(),
              ylim=c(0, 1),
              custom_colors=c("red", "blue"),
              custom_linetypes=c("solid", "dashed"),
              ci_draw_alpha=0.4,
              steps=TRUE,
              censoring_ind="lines",
              censoring_ind_width=0.1,
              censoring_ind_size=0.6)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, using many many things", fig=plt)
})
