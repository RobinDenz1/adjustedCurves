
library(riskRegression)

set.seed(42)
sim_dat <- sim_confounded_crisk(n=20)
sim_dat$group <- as.factor(sim_dat$group)

adj <- cif_aalen_johansen(data=sim_dat,
                          variable="group",
                          ev_time="time",
                          event="event",
                          conf_int=FALSE,
                          cause=1)

test_that("print.adjustedcif.method", {
  expect_snapshot_output(print(adj))
})

test_that("summary.adjustedcif.method", {
  expect_snapshot_output(summary(adj))
})
