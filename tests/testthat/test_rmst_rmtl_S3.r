library(survival)
library(ggplot2)
library(vdiffr)

set.seed(42)
sim_dat <- sim_confounded_surv(n=50, max_t=1.5)
sim_dat$group <- factor(sim_dat$group)

adj <- adjustedsurv(data=sim_dat,
                    variable="group",
                    ev_time="time",
                    event="event",
                    method="km",
                    conf_int=TRUE,
                    bootstrap=TRUE,
                    n_boot=10)

adj_rmst <- adjusted_rmst(adjsurv=adj, to=1.1, from=0, use_boot=TRUE)
adj_rmtl <- adjusted_rmtl(adj=adj, to=1.1, from=0, use_boot=TRUE)

## print methods

test_that("print.adjusted_rmst, default", {
  expect_snapshot_output(print(adj_rmst))
})

test_that("print.adjusted_rmst, digits", {
  expect_snapshot_output(print(adj_rmst, digits=2), NA)
})

test_that("print.adjusted_rmtl, default", {
  expect_snapshot_output(print(adj_rmtl), NA)
})

test_that("print.adjusted_rmtl, digits", {
  expect_snapshot_output(print(adj_rmtl, digits=2), NA)
})

## summary methods

test_that("summary.adjusted_rmst, default", {
  expect_snapshot_output(summary(adj_rmst), NA)
})

test_that("summary.adjusted_rmst, digits", {
  expect_snapshot_output(summary(adj_rmst, digits=2), NA)
})

test_that("summary.adjusted_rmtl, default", {
  expect_snapshot_output(print(adj_rmst), NA)
})

test_that("summary.adjusted_rmtl, digits", {
  expect_snapshot_output(print(adj_rmst, digits=2), NA)
})

## plot methods

test_that("plot.adjusted_rmst, default", {
  plt <- plot(adj_rmst)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot.adjusted_rmst, default", fig=plt)
})

test_that("plot.adjusted_rmst, conf_int", {
  plt <- plot(adj_rmst, conf_int=FALSE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot.adjusted_rmst, conf_int", fig=plt)
})

test_that("plot.adjusted_rmst, color", {
  plt <- plot(adj_rmst, color=FALSE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot.adjusted_rmst, color", fig=plt)
})

test_that("plot.adjusted_rmst, point_size", {
  plt <- plot(adj_rmst, point_size=4)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot.adjusted_rmst, point_size", fig=plt)
})

test_that("plot.adjusted_rmst, ci_size", {
  plt <- plot(adj_rmst, ci_size=2)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot.adjusted_rmst, ci_size", fig=plt)
})

test_that("plot.adjusted_rmst, ci_width", {
  plt <- plot(adj_rmst, ci_width=0.9)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot.adjusted_rmst, ci_width", fig=plt)
})

test_that("plot.adjusted_rmst, xlab + ylab + title", {
  plt <- plot(adj_rmst, xlab="X", ylab="Y", title="Title")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot.adjusted_rmst, xlab + ylab + title",
                              fig=plt)
})

test_that("plot.adjusted_rmst, ggplot theme", {
  plt <- plot(adj_rmst, gg_theme=ggplot2::theme_bw())
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot.adjusted_rmst, ggplot theme", fig=plt)
})

test_that("plot.adjusted_rmtl, default", {
  plt <- plot(adj_rmtl)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot.adjusted_rmtl, default", fig=plt)
})

test_that("plot.adjusted_rmtl, conf_int", {
  plt <- plot(adj_rmtl, conf_int=FALSE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot.adjusted_rmtl, conf_int", fig=plt)
})

test_that("plot.adjusted_rmtl, color", {
  plt <- plot(adj_rmtl, color=FALSE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot.adjusted_rmtl, color", fig=plt)
})

test_that("plot.adjusted_rmtl, point_size", {
  plt <- plot(adj_rmtl, point_size=4)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot.adjusted_rmtl, point_size", fig=plt)
})

test_that("plot.adjusted_rmtl, ci_size", {
  plt <- plot(adj_rmtl, ci_size=2)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot.adjusted_rmtl, ci_size", fig=plt)
})

test_that("plot.adjusted_rmtl, ci_width", {
  plt <- plot(adj_rmtl, ci_width=0.9)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot.adjusted_rmtl, ci_width", fig=plt)
})

test_that("plot.adjusted_rmtl, xlab + ylab + title", {
  plt <- plot(adj_rmtl, xlab="X", ylab="Y", title="Title")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot.adjusted_rmtl, xlab + ylab + title",
                              fig=plt)
})

test_that("plot.adjusted_rmtl, ggplot theme", {
  plt <- plot(adj_rmtl, gg_theme=ggplot2::theme_bw())
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot.adjusted_rmtl, ggplot theme", fig=plt)
})
