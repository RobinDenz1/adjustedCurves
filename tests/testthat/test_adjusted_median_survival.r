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

test_that("median surv, no boot", {
  expect_error(adjusted_median_survival(adj, use_boot=FALSE), NA)
})

test_that("median surv, with boot", {
  expect_error(adjusted_median_survival(adj, use_boot=TRUE), NA)
})

