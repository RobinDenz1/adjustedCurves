
suppressMessages(requireNamespace("survival"))

set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_crisk_n_50.Rds",
                               package="adjustedCurves"))
sim_dat$group <- as.factor(sim_dat$group)

## Just check if function throws any errors
test_that("2 treatments, no conf_int, no boot, no ...", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="aalen_johansen",
                     conf_int=FALSE,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

test_that("defaults, specific times", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="aalen_johansen",
                     bootstrap=TRUE,
                     n_boot=2,
                     conf_int=TRUE,
                     times=c(0, 0.2, 0.4),
                     cause=1)
  p1 <- plot(adj)
  p2 <- plot(adj, conf_int=TRUE, use_boot=TRUE)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

test_that("2 treatments, with conf_int, no boot, no ...", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="aalen_johansen",
                     conf_int=TRUE,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

test_that("2 treatments, no conf_int, with boot, no ...", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="aalen_johansen",
                     conf_int=FALSE,
                     bootstrap=TRUE,
                     n_boot=2,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

test_that("2 treatments, with conf_int, with boot, no ...", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="aalen_johansen",
                     conf_int=TRUE,
                     bootstrap=TRUE,
                     n_boot=2,
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


test_that("> 2 treatments, no conf_int, no boot, no ...", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="aalen_johansen",
                     conf_int=FALSE,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

test_that("> 2 treatments, with conf_int, no boot, no ...", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="aalen_johansen",
                     conf_int=TRUE,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

test_that("> 2 treatments, no conf_int, with boot, no ...", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="aalen_johansen",
                     conf_int=FALSE,
                     bootstrap=TRUE,
                     n_boot=2,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

test_that("> 2 treatments, with conf_int, with boot, no ...", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="aalen_johansen",
                     conf_int=TRUE,
                     bootstrap=TRUE,
                     n_boot=2,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})
