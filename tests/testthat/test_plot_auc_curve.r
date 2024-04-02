
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

## RMST
test_that("plot, rmst no arguments", {
  plt <- plot_rmst_curve(adj)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, rmst no arguments", fig=plt)
})

test_that("plot, rmst conf_int", {
  plt <- plot_rmst_curve(adj, conf_int=TRUE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, rmst conf_int", fig=plt)
})

test_that("plot difference, rmst", {
  plt <- plot_rmst_curve(adj, contrast="diff")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot difference, rmst", fig=plt)
})

test_that("plot difference, rmst conf_int", {
  plt <- plot_rmst_curve(adj, contrast="diff", conf_int=TRUE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot difference, rmst conf_int", fig=plt)
})

test_that("plot ratio, rmst", {
  plt <- plot_rmst_curve(adj, contrast="ratio")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot ratio, rmst", fig=plt)
})

test_that("plot ratio, rmst conf_int", {
  plt <- plot_rmst_curve(adj, contrast="ratio", conf_int=TRUE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot ratio, rmst conf_int", fig=plt)
})

# technically still only RMST, but general stuff
test_that("plot, max_t", {
  plt <- plot_rmst_curve(adj, max_t=0.4)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, max_t", fig=plt)
})

test_that("plot, times", {
  plt <- plot_rmst_curve(adj, times=seq(0.1, 0.8, 0.05))
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, times", fig=plt)
})

test_that("plot, linetype", {
  plt <- plot_rmst_curve(adj, linetype=TRUE, color=FALSE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, linetype", fig=plt)
})

test_that("plot, facet", {
  plt <- plot_rmst_curve(adj, facet=TRUE, color=FALSE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, facet", fig=plt)
})

test_that("plot, custom_linetypes", {
  plt <- plot_rmst_curve(adj, linetype=TRUE, color=FALSE,
                         custom_linetypes=c("dotdash", "dashed"))
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, custom_linetypes", fig=plt)
})

test_that("plot, labs", {
  plt <- plot_rmst_curve(adj, xlab="X", ylab="Y", title="A", subtitle="B",
                         legend.title="Legend")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, labs", fig=plt)
})

test_that("plot, linetype", {
  plt <- plot_rmst_curve(adj, linetype=TRUE, color=FALSE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, linetype", fig=plt)
})

test_that("plot difference, max_t", {
  plt <- plot_rmst_curve(adj, max_t=0.4, contrast="diff")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot difference, max_t", fig=plt)
})

test_that("plot difference, times", {
  plt <- plot_rmst_curve(adj, times=seq(0.1, 0.8, 0.05), contrast="diff")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot difference, times", fig=plt)
})

test_that("plot ratio, color", {
  plt <- plot_rmst_curve(adj, contrast="ratio", color="blue")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot ratio, color", fig=plt)
})

test_that("plot ratio, linetype", {
  plt <- plot_rmst_curve(adj, contrast="ratio", linetype="dashed")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot ratio, linetype", fig=plt)
})

## RMTL
test_that("plot, rmtl no arguments", {
  plt <- plot_rmtl_curve(adj)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, rmtl no arguments", fig=plt)
})

test_that("plot, rmtl conf_int", {
  plt <- plot_rmtl_curve(adj, conf_int=TRUE, color=FALSE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, rmtl conf_int", fig=plt)
})

test_that("plot difference, rmtl no arguments", {
  plt <- plot_rmtl_curve(adj, contrast="diff")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot difference, rmtl no arguments", fig=plt)
})

test_that("plot difference, rmtl conf_int", {
  plt <- plot_rmtl_curve(adj, conf_int=TRUE, contrast="diff")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot difference, rmtl conf_int", fig=plt)
})

## error messages
test_that("wrong adj rmtl", {
  expect_error(plot_rmtl_curve("string"),
               paste0("'adj' must be either an adjustedsurv object created",
                      " using the adjustedsurv function or an adjustedcif",
                      " object created using the adjustedcif function."),
               fixed=TRUE)
})

test_that("wrong adjsurv rmst", {
  expect_error(plot_rmst_curve("string"),
               paste0("'adjsurv' must be an adjustedsurv object created ",
                      "using the adjustedsurv function."),
               fixed=TRUE)
})
