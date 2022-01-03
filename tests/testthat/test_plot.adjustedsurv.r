
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
              censoring_ind_shape=10,
              censoring_ind_size=5,
              censoring_ind_alpha=0.5)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, using censoring indicators (points)",
                              fig=plt)
})

test_that("plot, using censoring indicators (lines) + things", {
  plt <- plot(adj, iso_reg=TRUE,
              censoring_ind="lines",
              censoring_ind_width=0.1,
              censoring_ind_size=1,
              censoring_ind_alpha=0.5,
              color=FALSE,
              single_color="blue",
              single_linetype="dashed",
              ylim=c(-0.4, 1.2))
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger(paste0("plot, using cens ind ",
                                     "(lines) and things"),
                              fig=plt)
})

test_that("plot, using censoring indicators (points) + things", {
  plt <- plot(adj, iso_reg=TRUE,
              censoring_ind="points",
              censoring_ind_shape=10,
              censoring_ind_size=5,
              censoring_ind_alpha=0.5,
              color=FALSE,
              single_color="blue",
              single_linetype="dashed",
              ylim=c(-0.4, 1.2))
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger(paste0("plot, using cens ind ",
                                     "(points) and things"),
                              fig=plt)
})

test_that("plot, using median surv lines", {
  plt <- plot(adj, iso_reg=TRUE,
              median_surv_lines=TRUE,
              median_surv_size=1,
              median_surv_linetype="solid",
              median_surv_color="red")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, using median surv lines", fig=plt)
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

test_that("plot, using single_color", {
  plt <- plot(adj, color=FALSE, single_color="blue")
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, using single_color", fig=plt)
})

test_that("plot, using single_linetype", {
  plt <- plot(adj, linetype=FALSE, color=FALSE, single_linetype="blue")
  expect_s3_class(plt, "ggplot")
})

test_that("plot, using conf_int_alpha", {
  plt <- plot(adj, conf_int=TRUE, conf_int_alpha=0.1)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, using conf_int_alpha", fig=plt)
})

test_that("plot, using steps", {
  plt <- suppressWarnings(plot(adj, steps=FALSE))
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, using steps", fig=plt)
})

test_that("plot, using cif", {
  plt <- plot(adj, cif=TRUE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, using cif", fig=plt)
})

test_that("plot, using no colors ci", {
  plt <- plot(adj, conf_int=TRUE, color=FALSE)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, using no colors ci", fig=plt)
})

test_that("plot, using no colors ci with steps", {
  plt <- suppressWarnings(plot(adj, conf_int=TRUE, color=FALSE, steps=FALSE))
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, ci no color steps", fig=plt)
})

test_that("plot, using single colors ci with steps", {
  plt <- suppressWarnings(plot(adj, conf_int=TRUE, color=FALSE, steps=FALSE,
                               single_color="red"))
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, ci single color steps", fig=plt)
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
              ylim=c(-0.1, 1.1),
              custom_colors=c("red", "blue"),
              custom_linetypes=c("solid", "dashed"),
              conf_int_alpha=0.4,
              steps=TRUE,
              median_surv_lines=TRUE,
              median_surv_size=1.2,
              median_surv_linetype="solid",
              median_surv_color="red",
              censoring_ind="lines",
              censoring_ind_width=0.1,
              censoring_ind_size=0.6)
  expect_s3_class(plt, "ggplot")
  vdiffr::expect_doppelganger("plot, using many many things", fig=plt)
})

## warnings and errors

test_that("Isotonic Regression with missing values", {
  adj_err <- adj
  adj_err$adjsurv$surv[1] <- NA
  expect_error(plot(adj_err, iso_reg=TRUE),
               paste0("Isotonic Regression cannot be used when there are ",
                      "missing values in the final survival estimates."))
})

test_that("single_color overwriting color", {
  expect_warning(plot(adj, color=TRUE, single_color="red"),
                paste0("Argument 'color' will be overwritten by ",
                       "argument 'single_color'."))
})

test_that("single_linetype overwriting linetype", {
  expect_warning(plot(adj, linetype=TRUE, single_linetype="dashed"),
                 paste0("Argument 'linetype' will be overwritten by ",
                        "argument 'single_linetype'."))
})

test_that("undefined censoring_ind argument", {
  expect_error(plot(adj, censoring_ind="undefined"),
                 paste0("Argument 'censoring_ind' must be either 'none', ",
                        "'lines' or 'points'. See documentation."))
})

test_that("use_boot with no boot no ci", {
  adj_err <- adj
  adj_err$boot_adjsurv <- NULL
  expect_warning(plot(adj_err, use_boot=TRUE, conf_int=TRUE),
                 paste0("Cannot use bootstrapped estimates as they were not ",
                        "estimated. Need bootstrap=TRUE in ",
                        "adjustedsurv() call."),
                 fixed=TRUE)
})

test_that("use_boot would work, conf_int not", {
  adj_err <- adj
  adj_err$adjsurv$ci_lower <- NULL
  expect_warning(plot(adj_err, use_boot=FALSE, conf_int=TRUE),
                 paste0("Cannot draw confidence intervals. Either set ",
                        "'conf_int=TRUE' in adjustedsurv() call or ",
                        "use bootstrap estimates."),
                 fixed=TRUE)
})

test_that("cif and median_surv_lines", {
  expect_warning(plot(adj, cif=TRUE, median_surv_lines=TRUE),
                 paste0("Cannot draw median survival indicators when ",
                        "using cif=TRUE."))
})
