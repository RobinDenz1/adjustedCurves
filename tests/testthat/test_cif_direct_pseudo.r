
set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_crisk_n_50.Rds",
                               package="adjustedCurves"))
sim_dat$group <- as.factor(sim_dat$group)

outcome_vars <- c("x1", "x2", "x3", "x4")

## Just check if function throws any errors
test_that("2 treatments, no conf_int, no boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct_pseudo",
                     conf_int=FALSE,
                     outcome_vars=outcome_vars,
                     type_time="bs",
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

test_that("2 treatments, no conf_int, with boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct_pseudo",
                     conf_int=FALSE,
                     bootstrap=TRUE,
                     n_boot=2,
                     outcome_vars=outcome_vars,
                     type_time="bs",
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

test_that("2 treatments, no conf_int, no boot, with times", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct_pseudo",
                     conf_int=FALSE,
                     bootstrap=TRUE,
                     n_boot=2,
                     outcome_vars=outcome_vars,
                     times=c(0.3, 0.8),
                     type_time="factor",
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

test_that(
  "2 treatments, no conf_int, no boot, with times, type_time='factor'", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct_pseudo",
                     conf_int=FALSE,
                     bootstrap=FALSE,
                     n_boot=5,
                     outcome_vars=outcome_vars,
                     times=c(0.3, 0.8),
                     type_time="factor",
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

test_that(
  "2 treatments, no conf_int, no boot, with times, type_time='ns'", {
    adj <- adjustedcif(data=sim_dat,
                       variable="group",
                       ev_time="time",
                       event="event",
                       method="direct_pseudo",
                       conf_int=FALSE,
                       bootstrap=FALSE,
                       n_boot=5,
                       outcome_vars=outcome_vars,
                       times=c(0.3, 0.8, 0.81),
                       type_time="ns",
                       spline_df=2,
                       cause=1)
    expect_s3_class(adj, "adjustedcif")
    expect_true(is.numeric(adj$adjcif$cif))
    expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
  })

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_crisk_n_50.Rds",
                               package="adjustedCurves"))
sim_dat$group[sim_dat$group==1] <- sample(c(1, 2),
                                        size=nrow(sim_dat[sim_dat$group==1, ]),
                                        replace=TRUE)
sim_dat$group <- as.factor(sim_dat$group)
sim_dat$time <- round(sim_dat$time, 1)

test_that("> 2 treatments, no conf_int, no boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct_pseudo",
                     conf_int=FALSE,
                     outcome_vars=outcome_vars,
                     type_time="bs",
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

test_that("> 2 treatments, no conf_int, with boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct_pseudo",
                     conf_int=FALSE,
                     bootstrap=TRUE,
                     n_boot=2,
                     outcome_vars=outcome_vars,
                     type_time="bs",
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

test_that("> 2 treatments, no conf_int, no boot, with times", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct_pseudo",
                     conf_int=FALSE,
                     bootstrap=FALSE,
                     n_boot=2,
                     outcome_vars=outcome_vars,
                     times=c(0.3, 0.8, 1, 1.5),
                     type_time="factor",
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

test_that("> 2 treatments, no conf_int, no boot, with times, type_time='factor'"
          , {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct_pseudo",
                     conf_int=FALSE,
                     bootstrap=FALSE,
                     n_boot=5,
                     outcome_vars=outcome_vars,
                     times=c(0.3, 0.8, 1, 1.5),
                     type_time="factor",
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})
