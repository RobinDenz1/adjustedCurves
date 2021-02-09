library(adjustedCurves)

set.seed(42)

sim_dat <- adjustedCurves::sim_confounded_surv(n=300)
sim_dat$group <- as.factor(sim_dat$group)

## Just check if function throws any errors
test_that("2 treatments, no sd, no boot, no ...", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="km",
                                            sd=F), NA)
})

test_that("2 treatments, with sd, no boot, no ...", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="km",
                                            sd=T), NA)
})

test_that("2 treatments, no sd, with boot, no ...", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="km",
                                            sd=F,
                                            bootstrap=T,
                                            n_boot=10), NA)
})

test_that("2 treatments, with sd, with boot, no ...", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="km",
                                            sd=T,
                                            bootstrap=T,
                                            n_boot=10), NA)
})

test_that("2 treatments, no sd, no boot, with ...", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="km",
                                            sd=F,
                                            error="greenwood"), NA)
})

test_that("2 treatments, with sd, with boot, with ...", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="km",
                                            sd=T,
                                            bootstrap=T,
                                            n_boot=10,
                                            error="greenwood"), NA)
})


sim_dat <- adjustedCurves::sim_confounded_surv(n=300)
sim_dat$group[sim_dat$group==1] <- sample(c(1, 2),
                                          size=nrow(sim_dat[sim_dat$group==1,]),
                                          replace=T)
sim_dat$group <- as.factor(sim_dat$group)


test_that("> 2 treatments, no sd, no boot, no ...", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="km",
                                            sd=F), NA)
})

test_that("> 2 treatments, with sd, no boot, no ...", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="km",
                                            sd=T), NA)
})

test_that("> 2 treatments, no sd, with boot, no ...", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="km",
                                            sd=F,
                                            bootstrap=T,
                                            n_boot=10), NA)
})

test_that("> 2 treatments, with sd, with boot, no ...", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="km",
                                            sd=T,
                                            bootstrap=T,
                                            n_boot=10), NA)
})

test_that("> 2 treatments, no sd, no boot, with ...", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="km",
                                            sd=F,
                                            error="greenwood"), NA)
})

test_that("> 2 treatments, with sd, with boot, with ...", {
  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
                                            variable="group",
                                            ev_time="time",
                                            event="event",
                                            method="km",
                                            sd=T,
                                            bootstrap=T,
                                            n_boot=10,
                                            error="greenwood"), NA)
})




