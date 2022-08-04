library(survival)
library(survtmle)

set.seed(35)

sim_dat <- sim_confounded_surv(n=11, max_t=1.5)
sim_dat$group <- as.factor(sim_dat$group)
sim_dat$time <- round(sim_dat$time * 10) + 1

# outcome model
adjust_vars <- c("x1", "x2", "x3", "x4", "x5", "x6")

test_that("2 treatments, no conf_int, no boot, using SL", {
  adj <- suppressWarnings(adjustedsurv(data=sim_dat,
                                       variable="group",
                                       ev_time="time",
                                       event="event",
                                       method="tmle",
                                       conf_int=FALSE,
                                       SL.trt=c("SL.glm"),
                                       SL.ftime=c("SL.glm"),
                                       SL.ctime=c("SL.glm"),
                                       adjust_vars=NULL))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), c("0", "1"))
})

test_that("2 treatments, no conf_int, with boot, using SL", {
  adj <- suppressWarnings(adjustedsurv(data=sim_dat,
                                       variable="group",
                                       ev_time="time",
                                       event="event",
                                       method="tmle",
                                       conf_int=FALSE,
                                       bootstrap=TRUE,
                                       n_boot=2,
                                       SL.trt=c("SL.glm"),
                                       SL.ftime=c("SL.glm"),
                                       SL.ctime=c("SL.glm"),
                                       adjust_vars=adjust_vars))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), c("0", "1"))
})

test_that("2 treatments, with conf_int, no boot, using SL", {
  adj <- suppressWarnings(adjustedsurv(data=sim_dat,
                                       variable="group",
                                       ev_time="time",
                                       event="event",
                                       method="tmle",
                                       conf_int=TRUE,
                                       bootstrap=FALSE,
                                       n_boot=2,
                                       SL.trt=c("SL.glm"),
                                       SL.ftime=c("SL.glm"),
                                       SL.ctime=NULL,
                                       adjust_vars=adjust_vars))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), c("0", "1"))
})

test_that("2 treatments, no conf_int, no boot, with times, using SL", {
  adj <- suppressWarnings(adjustedsurv(data=sim_dat,
                                       variable="group",
                                       ev_time="time",
                                       event="event",
                                       method="tmle",
                                       conf_int=FALSE,
                                       bootstrap=FALSE,
                                       SL.trt=c("SL.glm"),
                                       SL.ftime=c("SL.glm"),
                                       SL.ctime=NULL,
                                       adjust_vars=adjust_vars,
                                       times=c(1, 2, 3, 4)))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), c("0", "1"))
})

test_that("2 treatments, no conf_int, no boot, using glm", {
  adj <- suppressWarnings(adjustedsurv(data=sim_dat,
                                       variable="group",
                                       ev_time="time",
                                       event="event",
                                       method="tmle",
                                       conf_int=FALSE,
                                       glm.trt="x2 + x3 + x4",
                                       glm.ftime="x2 + x4 + x6",
                                       glm.ctime="x1 + x6",
                                       adjust_vars=adjust_vars))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), c("0", "1"))
})

test_that("2 treatments, no conf_int, with boot, using glm", {
  adj <- suppressWarnings(adjustedsurv(data=sim_dat,
                                       variable="group",
                                       ev_time="time",
                                       event="event",
                                       method="tmle",
                                       conf_int=FALSE,
                                       bootstrap=TRUE,
                                       n_boot=2,
                                       glm.trt="x2 + x3 + x4",
                                       glm.ftime="x2 + x4 + x6",
                                       glm.ctime="x1 + x6",
                                       adjust_vars=adjust_vars))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), c("0", "1"))
})

test_that("2 treatments, with conf_int, no boot, using glm", {
  adj <- suppressWarnings(adjustedsurv(data=sim_dat,
                                       variable="group",
                                       ev_time="time",
                                       event="event",
                                       method="tmle",
                                       conf_int=TRUE,
                                       bootstrap=FALSE,
                                       n_boot=2,
                                       glm.trt="x2 + x3 + x4",
                                       glm.ftime="x2 + x4 + x6",
                                       glm.ctime="x1 + x6",
                                       adjust_vars=adjust_vars))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), c("0", "1"))
})

test_that("2 treatments, no conf_int, no boot, with times, using glm", {
  adj <- suppressWarnings(adjustedsurv(data=sim_dat,
                                       variable="group",
                                       ev_time="time",
                                       event="event",
                                       method="tmle",
                                       conf_int=FALSE,
                                       bootstrap=FALSE,
                                       glm.trt="x2 + x3 + x4",
                                       glm.ftime="x2 + x4 + x6",
                                       glm.ctime="x1 + x6",
                                       adjust_vars=adjust_vars,
                                       times=c(1, 2, 3, 4)))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), c("0", "1"))
})

test_that("2 treatments, no conf_int, no boot, using SL, factor group", {
  sim_dat$group <- as.factor(sim_dat$group)
  adj <- suppressWarnings(adjustedsurv(data=sim_dat,
                                       variable="group",
                                       ev_time="time",
                                       event="event",
                                       method="tmle",
                                       conf_int=FALSE,
                                       SL.trt=c("SL.glm"),
                                       SL.ftime=c("SL.glm"),
                                       SL.ctime=c("SL.glm"),
                                       adjust_vars=NULL))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), c("0", "1"))
})
