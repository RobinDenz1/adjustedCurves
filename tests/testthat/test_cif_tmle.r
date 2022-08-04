library(survtmle)

set.seed(35)

sim_dat <- sim_confounded_surv(n=31, max_t=1.5)
sim_dat$group <- factor(sim_dat$group)
sim_dat$time <- round(sim_dat$time * 10) + 1
sim_dat$event[sim_dat$event==1] <- sample(c(1, 2), size=sum(sim_dat$event),
                                          replace=TRUE)
# outcome model
adjust_vars <- c("x1", "x2", "x3", "x4", "x5", "x6")

## Just check if function throws any errors
test_that("2 treatments, no conf_int, no boot, using SL", {
  adj <- suppressWarnings(adjustedcif(data=sim_dat,
                                      variable="group",
                                      ev_time="time",
                                      event="event",
                                      method="tmle",
                                      conf_int=FALSE,
                                      SL.trt=c("SL.glm"),
                                      SL.ftime=c("SL.glm"),
                                      SL.ctime=c("SL.glm"),
                                      adjust_vars=adjust_vars,
                                      cause=1))
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), c("0", "1"))
})

test_that("2 treatments, no conf_int, with boot, using SL", {
  adj <- suppressWarnings(adjustedcif(data=sim_dat,
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
                                      adjust_vars=adjust_vars,
                                      cause=1))
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), c("0", "1"))
})

test_that("2 treatments, with conf_int, no boot, using SL", {
  adj <- suppressWarnings(adjustedcif(data=sim_dat,
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
                                      adjust_vars=adjust_vars,
                                      times=c(1, 2, 3, 4, 5),
                                      cause=1))
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), c("0", "1"))
})

test_that("2 treatments, no conf_int, no boot, with times, using SL", {
  adj <- suppressWarnings(adjustedcif(data=sim_dat,
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
                                      cause=1))
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), c("0", "1"))
})

test_that("2 treatments, no conf_int, no boot, using glm", {
  adj <- suppressWarnings(adjustedcif(data=sim_dat,
                                      variable="group",
                                      ev_time="time",
                                      event="event",
                                      method="tmle",
                                      conf_int=FALSE,
                                      glm.trt="x2 + x3 + x4",
                                      glm.ftime="x2 + x4 + x6",
                                      glm.ctime="x1 + x6",
                                      adjust_vars=adjust_vars,
                                      cause=1))
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), c("0", "1"))
})

test_that("2 treatments, no conf_int, with boot, using glm", {
  adj <- suppressWarnings(adjustedcif(data=sim_dat,
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
                                      adjust_vars=adjust_vars,
                                      cause=1))
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), c("0", "1"))
})

test_that("2 treatments, with conf_int, no boot, using glm", {
  adj <- suppressWarnings(adjustedcif(data=sim_dat,
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
                                      adjust_vars=adjust_vars,
                                      times=c(1, 2, 3, 4, 5, 6),
                                      cause=1))
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), c("0", "1"))
})

test_that("2 treatments, no conf_int, no boot, with times, using glm", {
  adj <- suppressWarnings(adjustedcif(data=sim_dat,
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
                                      times=c(1, 2, 3, 4),
                                      cause=1))
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), c("0", "1"))
})
