library(ggplot2)

set.seed(42)

sim_dat <- sim_confounded_surv(n=20)
sim_dat$group <- as.factor(sim_dat$group)

adj <- adjustedsurv(data=sim_dat,
                    variable="group",
                    ev_time="time",
                    event="event",
                    method="km",
                    bootstrap=TRUE,
                    n_boot=20)

adj_test <- test_curve_equality(adjsurv=adj, from=0, to=0.5)

# 2 treatments

test_that("plot.curve_test, 2 treatments, type='curves'", {
  expect_error(plot(adj_test, type="curves"), NA)
})

test_that("plot.curve_test, 2 treatments, type='curves' + labs", {
  expect_error(plot(adj_test,
                    type="curves",
                    xlab="X",
                    ylab="Y",
                    title="Title"), NA)
})

test_that("plot.curve_test, 2 treatments, type='integral'", {
  expect_error(plot(adj_test, type="integral"), NA)
})

test_that("plot.curve_test, 2 treatments, type='integral' + labs", {
  expect_error(plot(adj_test,
                    type="integral",
                    xlab="X",
                    ylab="Y",
                    title="Title"), NA)
})

# > 2 treatments

sim_dat$group <- as.character(sim_dat$group)
sim_dat$group[sim_dat$group=="1"] <- sample(c(1, 2),
                                      size=nrow(sim_dat[sim_dat$group=="1",]),
                                      replace=TRUE)
sim_dat$group <- factor(sim_dat$group)

adj <- adjustedsurv(data=sim_dat,
                    variable="group",
                    ev_time="time",
                    event="event",
                    method="km",
                    bootstrap=TRUE,
                    n_boot=20)

adj_test <- test_curve_equality(adjsurv=adj, from=0, to=0.5)

test_that("plot.curve_test, > 2 treatments, type='curves'", {
  expect_error(plot(adj_test, type="curves"), NA)
})

test_that("plot.curve_test, > 2 treatments, type='curves' + labs", {
  expect_error(plot(adj_test,
                    type="curves",
                    xlab="X",
                    ylab="Y",
                    title="Title"), NA)
})

test_that("plot.curve_test, > 2 treatments, type='integral'", {
  expect_error(plot(adj_test, type="integral"), NA)
})

test_that("plot.curve_test, > 2 treatments, type='integral' + labs", {
  expect_error(plot(adj_test,
                    type="integral",
                    xlab="X",
                    ylab="Y",
                    title="Title"), NA)
})
