library(survival)

set.seed(42)

sim_dat <- sim_confounded_surv(n=100)
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
  expect_error(plot(adj), NA)
})

test_that("plot, with conf_int", {
  expect_error(plot(adj, draw_ci=TRUE), NA)
})

test_that("plot, using boot", {
  expect_error(plot(adj, draw_ci=TRUE, use_boot=TRUE), NA)
})

test_that("plot, using iso_reg", {
  expect_error(plot(adj, iso_reg=TRUE), NA)
})

test_that("plot, using censoring indicators", {
  expect_error(plot(adj, iso_reg=TRUE, censoring_ind=TRUE,
                    censoring_ind_width=0.1,
                    censoring_ind_size=1), NA)
})

test_that("plot, using median surv lines", {
  expect_error(plot(adj, iso_reg=TRUE,
                    median_surv_lines=TRUE,
                    median_surv_size=1,
                    median_surv_linetype="solid",
                    median_surv_color="red"), NA)
})

test_that("plot, using max_t", {
  expect_error(plot(adj, max_t=0.2), NA)
})

test_that("plot, using force_bounds", {
  expect_error(plot(adj, force_bounds=TRUE), NA)
})

test_that("plot, using linetype + facets + color", {
  expect_error(plot(adj, color=TRUE, linetype=TRUE, facet=TRUE), NA)
})

test_that("plot, using labs + title", {
  expect_error(plot(adj, xlab="X", ylab="Y", title="Title"), NA)
})

test_that("plot, using legend.title + legend.position", {
  expect_error(plot(adj, legend.title="A", legend.position="bottom"), NA)
})

test_that("plot, using ggplot theme", {
  expect_error(plot(adj, gg_theme=ggplot2::theme_bw()), NA)
})

test_that("plot, using ylim", {
  expect_error(plot(adj, ylim=c(0, 1)), NA)
})

test_that("plot, using custom_colors", {
  expect_error(plot(adj, custom_colors=c("red", "blue")), NA)
})

test_that("plot, using custom_linetypes", {
  expect_error(plot(adj, custom_linetypes=c("dashed", "solid"),
                    linetype=TRUE), NA)
})

test_that("plot, using ci_draw_alpha", {
  expect_error(plot(adj, draw_ci=TRUE, ci_draw_alpha=0.1), NA)
})

test_that("plot, using steps", {
  expect_error(plot(adj, steps=FALSE), NA)
})

test_that("plot, using cif", {
  expect_error(plot(adj, cif=TRUE), NA)
})

test_that("plot, using many many things", {
  expect_error(plot(adj,
                    draw_ci=TRUE,
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
                    median_surv_lines=TRUE,
                    median_surv_size=1.2,
                    median_surv_linetype="solid",
                    median_surv_color="red",
                    censoring_ind=TRUE,
                    censoring_ind_width=0.1,
                    censoring_ind_size=0.6), NA)
})
