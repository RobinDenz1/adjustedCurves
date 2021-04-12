library(survival)
library(survtmle)

set.seed(35)

sim_dat <- adjustedCurves::sim_confounded_surv(n=200, max_t=1.5)
sim_dat$time <- round(sim_dat$time * 10) + 1

# outcome model
adjust_vars <- c("x1", "x2", "x3", "x4", "x5", "x6")

## Just check if function throws any errors
test_that("2 treatments, no conf_int, no boot", {
  expect_error(suppressWarnings(adjustedsurv(data=sim_dat,
                                             variable="group",
                                             ev_time="time",
                                             event="event",
                                             method="ostmle",
                                             conf_int=F,
                                             SL.trt=c("SL.glm"),
                                             SL.ftime=c("SL.glm"),
                                             SL.ctime=c("SL.glm"),
                                             adjust_vars=adjust_vars)), NA)
})

test_that("2 treatments, no conf_int, with boot", {
  expect_error(suppressWarnings(adjustedsurv(data=sim_dat,
                                             variable="group",
                                             ev_time="time",
                                             event="event",
                                             method="ostmle",
                                             conf_int=F,
                                             bootstrap=T,
                                             n_boot=2,
                                             SL.trt=c("SL.glm"),
                                             SL.ftime=c("SL.glm"),
                                             SL.ctime=c("SL.glm"),
                                             adjust_vars=adjust_vars)), NA)
})

test_that("2 treatments, with conf_int, no boot", {
  expect_error(suppressWarnings(adjustedsurv(data=sim_dat,
                                             variable="group",
                                             ev_time="time",
                                             event="event",
                                             method="tmle",
                                             conf_int=T,
                                             bootstrap=F,
                                             n_boot=2,
                                             SL.trt=c("SL.glm"),
                                             SL.ftime=c("SL.glm"),
                                             SL.ctime=NULL,
                                             adjust_vars=adjust_vars)), NA)
})

test_that("2 treatments, no conf_int, no boot, with times", {
  expect_error(suppressWarnings(adjustedsurv(data=sim_dat,
                                             variable="group",
                                             ev_time="time",
                                             event="event",
                                             method="tmle",
                                             conf_int=F,
                                             bootstrap=F,
                                             SL.trt=c("SL.glm"),
                                             SL.ftime=c("SL.glm"),
                                             SL.ctime=NULL,
                                             adjust_vars=adjust_vars,
                                             times=c(1, 2, 3, 4))), NA)
})

test_that("2 treatments, no conf_int, no boot, with times, epsilon", {
  expect_error(suppressWarnings(adjustedsurv(data=sim_dat,
                                             variable="group",
                                             ev_time="time",
                                             event="event",
                                             method="tmle",
                                             conf_int=F,
                                             bootstrap=F,
                                             SL.trt=c("SL.glm"),
                                             SL.ftime=c("SL.glm"),
                                             SL.ctime=NULL,
                                             adjust_vars=adjust_vars,
                                             times=c(1, 2, 3, 4),
                                             epsilon=1.1)), NA)
})

test_that("2 treatments, no conf_int, no boot, with times, n_iter", {
  expect_error(suppressWarnings(adjustedsurv(data=sim_dat,
                                             variable="group",
                                             ev_time="time",
                                             event="event",
                                             method="tmle",
                                             conf_int=F,
                                             bootstrap=F,
                                             SL.trt=c("SL.glm"),
                                             SL.ftime=c("SL.glm"),
                                             SL.ctime=NULL,
                                             adjust_vars=adjust_vars,
                                             times=c(1, 2, 3, 4),
                                             max_num_iteration=50)), NA)
})

test_that("2 treatments, no conf_int, no boot, with times, moss_method", {
  expect_error(suppressWarnings(adjustedsurv(data=sim_dat,
                                             variable="group",
                                             ev_time="time",
                                             event="event",
                                             method="tmle",
                                             conf_int=F,
                                             bootstrap=F,
                                             SL.trt=c("SL.glm"),
                                             SL.ftime=c("SL.glm"),
                                             SL.ctime=NULL,
                                             adjust_vars=adjust_vars,
                                             times=c(1, 2, 3, 4),
                                             psi_moss_method="l1")), NA)
})

test_that("2 treatments, no conf_int, no boot, with times, gtol", {
  expect_error(suppressWarnings(adjustedsurv(data=sim_dat,
                                             variable="group",
                                             ev_time="time",
                                             event="event",
                                             method="tmle",
                                             conf_int=F,
                                             bootstrap=F,
                                             SL.trt=c("SL.glm"),
                                             SL.ftime=c("SL.glm"),
                                             SL.ctime=NULL,
                                             adjust_vars=adjust_vars,
                                             times=c(1, 2, 3, 4),
                                             gtol=0.1)), NA)
})

