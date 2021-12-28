
set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_crisk_n_100.Rds",
                               package="adjustedCurves"))
sim_dat$group <- as.factor(sim_dat$group)

mod <- glm(group ~ x1 + x2 + x3 + x4 + x5 + x6, data=sim_dat,
           family="binomial")

## Just check if function throws any errors
test_that("2 treatments, no conf_int, no boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="iptw",
                     conf_int=FALSE,
                     treatment_model=mod,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

test_that("2 treatments, with conf_int, no boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="iptw",
                     conf_int=TRUE,
                     treatment_model=mod,
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
                     method="iptw",
                     conf_int=FALSE,
                     bootstrap=TRUE,
                     n_boot=2,
                     treatment_model=mod,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

test_that("2 treatments, with conf_int, with boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="iptw",
                     conf_int=TRUE,
                     bootstrap=TRUE,
                     n_boot=2,
                     treatment_model=mod,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adjcif$cif))
  expect_equal(levels(adj$adjcif$group), levels(sim_dat$group))
})

#sim_dat <- sim_confounded_crisk(n=150)
#sim_dat$group[sim_dat$group==1] <- sample(c(1, 2),
#                                        size=nrow(sim_dat[sim_dat$group==1, ]),
#                                        replace=TRUE)
#sim_dat$group <- as.factor(sim_dat$group)
#
#mod <- nnet::multinom(group ~ x1 + x2 + x3 + x4 + x5 + x6, data=sim_dat)

# DONT RUN:
# - There is a known issue with the ate function of the riskRegression
#   R-Package when using a multinom model object. I have reported this
#   error to the developers and they have fixed it immediatly.
#   The fix however has not yet been uploaded to CRAN. Using the github
#   release version it is possible to run this code, but not with the
#   latest CRAN version.
#test_that("> 2 treatments, no conf_int, no boot, no ...", {
#  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
#                                           variable="group",
#                                           ev_time="time",
#                                           event="event",
#                                           method="iptw",
#                                           conf_int=FALSE,
#                                           treatment_model=mod,
#                                           cause=1), NA)
#})
#
#test_that("> 2 treatments, with conf_int, no boot, no ...", {
#  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
#                                           variable="group",
#                                           ev_time="time",
#                                           event="event",
#                                           method="iptw",
#                                           conf_int=TRUE,
#                                           treatment_model=mod,
#                                           cause=1), NA)
#})
#
#test_that("> 2 treatments, no conf_int, with boot, no ...", {
#  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
#                                           variable="group",
#                                           ev_time="time",
#                                           event="event",
#                                           method="iptw",
#                                           conf_int=FALSE,
#                                           bootstrap=TRUE,
#                                           n_boot=2,
#                                           treatment_model=mod,
#                                           cause=1), NA)
#})
#
#test_that("> 2 treatments, with conf_int, with boot, no ...", {
#  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
#                                           variable="group",
#                                           ev_time="time",
#                                           event="event",
#                                           method="iptw",
#                                           conf_int=TRUE,
#                                           bootstrap=TRUE,
#                                           n_boot=2,
#                                           treatment_model=mod,
#                                           cause=1), NA)
#})
