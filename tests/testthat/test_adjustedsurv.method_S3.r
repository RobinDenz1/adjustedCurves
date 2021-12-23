
set.seed(42)
sim_dat <- sim_confounded_surv(n=20)
sim_dat$group <- as.factor(sim_dat$group)

adj <- surv_km(data=sim_dat,
               variable="group",
               ev_time="time",
               event="event",
               conf_int=FALSE)

test_that("print.adjustedsurv.method", {
  expect_snapshot_output(print(adj))
})

test_that("summary.adjustedsurv.method", {
  expect_snapshot_output(summary(adj))
})
