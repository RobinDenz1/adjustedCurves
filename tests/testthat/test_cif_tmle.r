library(survival)
library(survtmle)

set.seed(35)

sim_dat <- adjustedCurves::sim_confounded_surv(n=100, max_t=1.5)
sim_dat$time <- round(sim_dat$time * 10) + 1
sim_dat$event[sim_dat$event==1] <- sample(c(1, 2), size=sum(sim_dat$event),
                                          replace=T)

# outcome model
adjust_vars <- c("x1", "x2", "x3", "x4", "x5", "x6")

## Just check if function throws any errors
test_that("2 treatments, no conf_int, no boot, using SL", {
  expect_error(suppressWarnings(adjustedcif(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="tmle",
                                            conf_int=F,
                                            SL.trt=c("SL.glm"),
                                            SL.ftime=c("SL.glm"),
                                            SL.ctime=c("SL.glm"),
                                            adjust_vars=adjust_vars,
                                            cause=1)), NA)
})

test_that("2 treatments, no conf_int, with boot, using SL", {
  expect_error(suppressWarnings(adjustedcif(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="tmle",
                                            conf_int=F,
                                            bootstrap=T,
                                            n_boot=2,
                                            SL.trt=c("SL.glm"),
                                            SL.ftime=c("SL.glm"),
                                            SL.ctime=c("SL.glm"),
                                            adjust_vars=adjust_vars,
                                            cause=1)), NA)
})

test_that("2 treatments, with conf_int, no boot, using SL", {
  expect_error(suppressWarnings(adjustedcif(data=sim_dat,
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
                                            adjust_vars=adjust_vars,
                                            cause=1)), NA)
})

test_that("2 treatments, no conf_int, no boot, with times, using SL", {
  expect_error(suppressWarnings(adjustedcif(data=sim_dat,
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
                                            cause=1)), NA)
})

test_that("2 treatments, no conf_int, no boot, using glm", {
  expect_error(suppressWarnings(adjustedcif(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="tmle",
                                            conf_int=F,
                                            glm.trt="x2 + x3 + x4",
                                            glm.ftime="x2 + x4 + x6",
                                            glm.ctime="x1 + x6",
                                            adjust_vars=adjust_vars,
                                            cause=1)), NA)
})

test_that("2 treatments, no conf_int, with boot, using glm", {
  expect_error(suppressWarnings(adjustedcif(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="tmle",
                                            conf_int=F,
                                            bootstrap=T,
                                            n_boot=2,
                                            glm.trt="x2 + x3 + x4",
                                            glm.ftime="x2 + x4 + x6",
                                            glm.ctime="x1 + x6",
                                            adjust_vars=adjust_vars,
                                            cause=1)), NA)
})

test_that("2 treatments, with conf_int, no boot, using glm", {
  expect_error(suppressWarnings(adjustedcif(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="tmle",
                                            conf_int=T,
                                            bootstrap=F,
                                            n_boot=2,
                                            glm.trt="x2 + x3 + x4",
                                            glm.ftime="x2 + x4 + x6",
                                            glm.ctime="x1 + x6",
                                            adjust_vars=adjust_vars,
                                            cause=1)), NA)
})

test_that("2 treatments, no conf_int, no boot, with times, using glm", {
  expect_error(suppressWarnings(adjustedcif(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="tmle",
                                            conf_int=F,
                                            bootstrap=F,
                                            glm.trt="x2 + x3 + x4",
                                            glm.ftime="x2 + x4 + x6",
                                            glm.ctime="x1 + x6",
                                            adjust_vars=adjust_vars,
                                            times=c(1, 2, 3, 4),
                                            cause=1)), NA)
})
