
set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_surv_n_50.Rds",
                               package="adjustedCurves"))
sim_dat$group <- as.factor(sim_dat$group)

test_that("2 treatments, no conf_int, no boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="km",
                      conf_int=FALSE)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("2 treatments, with conf_int, no boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="km",
                      conf_int=TRUE)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("2 treatments, no conf_int, with boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="km",
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=2)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("2 treatments, with conf_int, with boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="km",
                      conf_int=TRUE,
                      bootstrap=TRUE,
                      n_boot=2)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_surv_n_50.Rds",
                               package="adjustedCurves"))
sim_dat$group[sim_dat$group==1] <- sample(c(1, 2),
                                        size=nrow(sim_dat[sim_dat$group==1, ]),
                                        replace=TRUE)
sim_dat$group <- as.factor(sim_dat$group)


test_that("> 2 treatments, no conf_int, no boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="km",
                      conf_int=FALSE)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("> 2 treatments, with conf_int, no boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="km",
                      conf_int=TRUE)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("> 2 treatments, no conf_int, with boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="km",
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=2)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("> 2 treatments, with conf_int, with boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="km",
                      conf_int=TRUE,
                      bootstrap=TRUE,
                      n_boot=2)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

# very little data
sim_dat <- sim_confounded_surv(n=25)
sim_dat$group <- as.factor(sim_dat$group)

test_that("little data, checking if add_zeros works", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="km",
                      conf_int=TRUE,
                      bootstrap=TRUE,
                      n_boot=2)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
  p <- suppressWarnings(plot(adj, median_surv_lines=TRUE))
  expect_s3_class(p, "ggplot")
  p <- suppressWarnings(plot(adj, use_boot=TRUE, conf_int=TRUE))
  expect_s3_class(p, "ggplot")
  p <- suppressWarnings(plot(adj, conf_int=TRUE))
  expect_s3_class(p, "ggplot")
})
