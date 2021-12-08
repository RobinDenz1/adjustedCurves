library(survival)
library(ggplot2)

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
  expect_error(print(adj_rmst), NA)
})

test_that("print.adjusted_rmst, digits", {
  expect_error(print(adj_rmst, digits=2), NA)
})

test_that("print.adjusted_rmtl, default", {
  expect_error(print(adj_rmtl), NA)
})

test_that("print.adjusted_rmtl, digits", {
  expect_error(print(adj_rmtl, digits=2), NA)
})

## summary methods

test_that("summary.adjusted_rmst, default", {
  expect_error(summary(adj_rmst), NA)
})

test_that("summary.adjusted_rmst, digits", {
  expect_error(summary(adj_rmst, digits=2), NA)
})

test_that("summary.adjusted_rmtl, default", {
  expect_error(print(adj_rmst), NA)
})

test_that("summary.adjusted_rmtl, digits", {
  expect_error(print(adj_rmst, digits=2), NA)
})

## plot methods

test_that("plot.adjusted_rmst, default", {
  expect_error(plot(adj_rmst), NA)
})

test_that("plot.adjusted_rmst, draw_ci", {
  expect_error(plot(adj_rmst, draw_ci=FALSE), NA)
})

test_that("plot.adjusted_rmst, color", {
  expect_error(plot(adj_rmst, color=FALSE), NA)
})

test_that("plot.adjusted_rmst, point_size", {
  expect_error(plot(adj_rmst, point_size=4), NA)
})

test_that("plot.adjusted_rmst, ci_size", {
  expect_error(plot(adj_rmst, ci_size=2), NA)
})

test_that("plot.adjusted_rmst, ci_width", {
  expect_error(plot(adj_rmst, ci_width=0.9), NA)
})

test_that("plot.adjusted_rmst, xlab + ylab + title", {
  expect_error(plot(adj_rmst, xlab="X", ylab="Y", title="Title"), NA)
})

test_that("plot.adjusted_rmst, ggplot theme", {
  expect_error(plot(adj_rmst, gg_theme=ggplot2::theme_bw()), NA)
})

test_that("plot.adjusted_rmtl, default", {
  expect_error(plot(adj_rmtl), NA)
})

test_that("plot.adjusted_rmtl, draw_ci", {
  expect_error(plot(adj_rmtl, draw_ci=FALSE), NA)
})

test_that("plot.adjusted_rmtl, color", {
  expect_error(plot(adj_rmtl, color=FALSE), NA)
})

test_that("plot.adjusted_rmtl, point_size", {
  expect_error(plot(adj_rmtl, point_size=4), NA)
})

test_that("plot.adjusted_rmtl, ci_size", {
  expect_error(plot(adj_rmtl, ci_size=2), NA)
})

test_that("plot.adjusted_rmtl, ci_width", {
  expect_error(plot(adj_rmtl, ci_width=0.9), NA)
})

test_that("plot.adjusted_rmtl, xlab + ylab + title", {
  expect_error(plot(adj_rmtl, xlab="X", ylab="Y", title="Title"), NA)
})

test_that("plot.adjusted_rmtl, ggplot theme", {
  expect_error(plot(adj_rmtl, gg_theme=ggplot2::theme_bw()), NA)
})

