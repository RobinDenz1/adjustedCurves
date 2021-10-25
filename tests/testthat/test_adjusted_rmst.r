library(survival)

set.seed(42)
sim_dat <- sim_confounded_surv(n=50, max_t=1.5)
sim_dat$group <- factor(sim_dat$group)

adj <- adjustedsurv(data=sim_dat,
                    variable="group",
                    ev_time="time",
                    event="event",
                    method="km",
                    conf_int=TRUE,
                    bootstrap=TRUE,
                    n_boot=10)

test_that("rmst, no boot", {
  expect_error(adjusted_rmst(adj, to=1.1), NA)
})

test_that("rmst, with boot", {
  expect_error(adjusted_rmst(adj, to=1.1, use_boot=TRUE), NA)
})

test_that("rmst, no boot, using from", {
  expect_error(adjusted_rmst(adj, to=1.1, from=0.3), NA)
})

test_that("rmst, with boot, using from", {
  expect_error(adjusted_rmst(adj, to=1.1, from=0.3, use_boot=TRUE), NA)
})
