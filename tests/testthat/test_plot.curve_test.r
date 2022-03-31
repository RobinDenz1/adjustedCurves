
set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_surv_n_20.Rds",
                               package="adjustedCurves"))
sim_dat$group <- as.factor(sim_dat$group)

adj <- adjustedsurv(data=sim_dat,
                    variable="group",
                    ev_time="time",
                    event="event",
                    method="km",
                    bootstrap=TRUE,
                    n_boot=20)

adj_test <- adjusted_curve_test(adj=adj, from=0, to=0.5)

# 2 treatments
# NOTE: I am not using snapshots + expect_doppelganger here
#       because it did not work very well with these kinds of plots

test_that("plot.curve_test, 2 treatments, type='curves'", {
  plt <- plot(adj_test, type="curves")
  expect_s3_class(plt, "ggplot")
})

test_that("plot.curve_test, 2 treatments, type='curves' + labs", {
  plt <- plot(adj_test,
              type="curves",
              xlab="X",
              ylab="Y",
              title="Title")
  expect_s3_class(plt, "ggplot")
})

test_that("plot.curve_test, 2 treatments, type='integral'", {
  plt <- plot(adj_test, type="integral")
  expect_s3_class(plt, "ggplot")
})

test_that("plot.curve_test, 2 treatments, type='integral' + labs", {
  plt <- plot(adj_test,
              type="integral",
              xlab="X",
              ylab="Y",
              title="Title")
  expect_s3_class(plt, "ggplot")
})

adj_test <- adjusted_curve_test(adj=adj, from=0, to=0.5, interpolation="linear")

test_that("plot.curve_test, 2 treatments, type='curves', linear", {
  plt <- plot(adj_test, type="curves")
  expect_s3_class(plt, "ggplot")
})

# > 2 treatments

sim_dat$group <- as.character(sim_dat$group)
sim_dat$group[sim_dat$group=="1"] <- sample(c(1, 2),
                                      size=nrow(sim_dat[sim_dat$group=="1", ]),
                                      replace=TRUE)
sim_dat$group <- factor(sim_dat$group)

adj <- adjustedsurv(data=sim_dat,
                    variable="group",
                    ev_time="time",
                    event="event",
                    method="km",
                    bootstrap=TRUE,
                    n_boot=20)

adj_test <- adjusted_curve_test(adj=adj, from=0, to=0.5)

test_that("plot.curve_test, > 2 treatments, type='curves'", {
  plt <- plot(adj_test, type="curves")
  expect_s3_class(plt, "ggplot")
})

test_that("plot.curve_test, > 2 treatments, type='curves' + labs", {
  plt <- plot(adj_test,
              type="curves",
              xlab="X",
              ylab="Y",
              title="Title")
  expect_s3_class(plt, "ggplot")
})

test_that("plot.curve_test, > 2 treatments, type='integral'", {
  plt <- plot(adj_test, type="integral")
  expect_s3_class(plt, "ggplot")
})

test_that("plot.curve_test, > 2 treatments, type='integral' + labs", {
  plt <- plot(adj_test,
              type="integral",
              xlab="X",
              ylab="Y",
              title="Title")
  expect_s3_class(plt, "ggplot")
})
