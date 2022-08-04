library(survival)
library(survtmle)

set.seed(35)

sim_dat <- sim_confounded_surv(n=35, max_t=1.5)
sim_dat$group <- factor(sim_dat$group)
sim_dat$time <- round(sim_dat$time * 10) + 1

# outcome model
adjust_vars <- c("x1", "x2", "x3", "x4", "x5", "x6")

test_that("2 treatments, no conf_int, no boot", {
  adj <- suppressWarnings(adjustedsurv(data=sim_dat,
                                       variable="group",
                                       ev_time="time",
                                       event="event",
                                       method="ostmle",
                                       conf_int=FALSE,
                                       SL.trt=c("SL.glm"),
                                       SL.ftime=c("SL.glm"),
                                       SL.ctime=c("SL.glm"),
                                       adjust_vars=NULL))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), c("0", "1"))
})

test_that("2 treatments, no conf_int, with boot", {
  adj <- suppressWarnings(adjustedsurv(data=sim_dat,
                                       variable="group",
                                       ev_time="time",
                                       event="event",
                                       method="ostmle",
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

test_that("2 treatments, with conf_int, no boot", {
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

test_that("2 treatments, no conf_int, no boot, with times", {
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

test_that("2 treatments, no conf_int, no boot, with times, epsilon", {
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
                                       times=c(1, 2, 3, 4),
                                       epsilon=1.1))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), c("0", "1"))
})

test_that("2 treatments, no conf_int, no boot, with times, n_iter", {
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
                                       times=c(1, 2, 3, 4),
                                       max_num_iteration=50))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), c("0", "1"))
})

test_that("2 treatments, no conf_int, no boot, with times, moss_method", {
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
                                       times=c(1, 2, 3, 4),
                                       psi_moss_method="l1"))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), c("0", "1"))
})

test_that("2 treatments, no conf_int, no boot, with times, gtol", {
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
                                       times=c(1, 2, 3, 4),
                                       gtol=0.1))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), c("0", "1"))
})

#test_that("2 treatments, with conf_int, no boot", {
#  adj <- suppressWarnings(adjustedsurv(data=sim_dat,
#                                       variable="group",
#                                       ev_time="time",
#                                       event="event",
#                                       method="ostmle",
#                                       conf_int=TRUE,
#                                       SL.trt=c("SL.glm"),
#                                       SL.ftime=c("SL.glm"),
#                                       SL.ctime=c("SL.glm"),
#                                       adjust_vars=NULL))
#  expect_s3_class(adj, "adjustedsurv")
#  expect_true(is.numeric(adj$adjsurv$surv))
#  expect_equal(levels(adj$adjsurv$group), c("0", "1"))
#})

test_that("2 treatments, psi_moss_method='glm'", {
  adj <- suppressWarnings(adjustedsurv(data=sim_dat,
                                       variable="group",
                                       ev_time="time",
                                       event="event",
                                       method="ostmle",
                                       conf_int=FALSE,
                                       SL.trt=c("SL.glm"),
                                       SL.ftime=c("SL.glm"),
                                       SL.ctime=c("SL.glm"),
                                       adjust_vars=NULL,
                                       psi_moss_method="glm"))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), c("0", "1"))
})

test_that("2 treatments, psi_moss_method='l1'", {
  adj <- suppressWarnings(adjustedsurv(data=sim_dat,
                                       variable="group",
                                       ev_time="time",
                                       event="event",
                                       method="ostmle",
                                       conf_int=FALSE,
                                       SL.trt=c("SL.glm"),
                                       SL.ftime=c("SL.glm"),
                                       SL.ctime=c("SL.glm"),
                                       adjust_vars=NULL,
                                       psi_moss_method="l1"))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), c("0", "1"))
})
