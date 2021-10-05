
set.seed(42)

sim_dat <- adjustedCurves::sim_confounded_surv(n=50)
sim_dat$event[sim_dat$event==1] <- sample(c(1, 2), size=sum(sim_dat$event),
                                          replace=T)
sim_dat$group <- as.factor(sim_dat$group)

## Just check if function throws any errors
test_that("2 treatments, no conf_int, no boot, no ...", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aalen_johansen",
                                           conf_int=F,
                                           cause=1), NA)
})

test_that("2 treatments, with conf_int, no boot, no ...", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aalen_johansen",
                                           conf_int=T,
                                           cause=1), NA)
})

test_that("2 treatments, no conf_int, with boot, no ...", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aalen_johansen",
                                           conf_int=F,
                                           bootstrap=T,
                                           n_boot=2,
                                           cause=1), NA)
})

test_that("2 treatments, with conf_int, with boot, no ...", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aalen_johansen",
                                           conf_int=T,
                                           bootstrap=T,
                                           n_boot=2,
                                           cause=1), NA)
})

sim_dat <- adjustedCurves::sim_confounded_surv(n=50)
sim_dat$group[sim_dat$group==1] <- sample(c(1, 2),
                                          size=nrow(sim_dat[sim_dat$group==1,]),
                                          replace=T)
sim_dat$event[sim_dat$event==1] <- sample(c(1, 2), size=sum(sim_dat$event),
                                          replace=T)
sim_dat$group <- as.factor(sim_dat$group)


test_that("> 2 treatments, no conf_int, no boot, no ...", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aalen_johansen",
                                           conf_int=F,
                                           cause=1), NA)
})

test_that("> 2 treatments, with conf_int, no boot, no ...", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aalen_johansen",
                                           conf_int=T,
                                           cause=1), NA)
})

test_that("> 2 treatments, no conf_int, with boot, no ...", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aalen_johansen",
                                           conf_int=F,
                                           bootstrap=T,
                                           n_boot=2,
                                           cause=1), NA)
})

test_that("> 2 treatments, with conf_int, with boot, no ...", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="aalen_johansen",
                                           conf_int=T,
                                           bootstrap=T,
                                           n_boot=2,
                                           cause=1), NA)
})
