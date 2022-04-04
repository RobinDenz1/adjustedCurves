
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
                    bootstrap=TRUE,
                    n_boot=10)

adj_rmst <- adjusted_rmst(adjsurv=adj, to=1.1, from=0, conf_int=TRUE)
adj_rmtl <- adjusted_rmtl(adj=adj, to=1.1, from=0, conf_int=TRUE)

adj_rmst_no_ci <- adjusted_rmst(adjsurv=adj, to=1.1, from=0, conf_int=FALSE)
adj_rmtl_no_ci <- adjusted_rmtl(adj=adj, to=1.1, from=0, conf_int=FALSE)

## print methods

test_that("print.adjusted_rmst no boot, default", {
  expect_snapshot_output(print(adj_rmst_no_ci))
})

test_that("print.adjusted_rmst no boot, digits", {
  expect_snapshot_output(print(adj_rmst_no_ci, digits=2))
})

test_that("print.adjusted_rmtl no boot, default", {
  expect_snapshot_output(print(adj_rmtl_no_ci))
})

test_that("print.adjusted_rmtl no boot, digits", {
  expect_snapshot_output(print(adj_rmtl_no_ci, digits=2))
})

test_that("print.adjusted_rmst with boot, default", {
  expect_snapshot_output(print(adj_rmst))
})

test_that("print.adjusted_rmst with boot, digits", {
  expect_snapshot_output(print(adj_rmst, digits=2))
})

test_that("print.adjusted_rmtl with boot, default", {
  expect_snapshot_output(print(adj_rmtl))
})

test_that("print.adjusted_rmtl with boot, digits", {
  expect_snapshot_output(print(adj_rmtl, digits=2))
})

## summary methods

test_that("summary.adjusted_rmst no boot, default", {
  expect_snapshot_output(summary(adj_rmst_no_ci))
})

test_that("summary.adjusted_rmst no boot, digits", {
  expect_snapshot_output(summary(adj_rmst_no_ci, digits=2))
})

test_that("summary.adjusted_rmtl no boot, default", {
  expect_snapshot_output(summary(adj_rmtl_no_ci))
})

test_that("summary.adjusted_rmtl no boot, digits", {
  expect_snapshot_output(summary(adj_rmtl_no_ci, digits=2))
})

test_that("summary.adjusted_rmst with boot, default", {
  expect_snapshot_output(summary(adj_rmst))
})

test_that("summary.adjusted_rmst with boot, digits", {
  expect_snapshot_output(summary(adj_rmst, digits=2))
})

test_that("summary.adjusted_rmtl with boot, default", {
  expect_snapshot_output(summary(adj_rmtl))
})

test_that("summary.adjusted_rmtl with boot, digits", {
  expect_snapshot_output(summary(adj_rmtl, digits=2))
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
