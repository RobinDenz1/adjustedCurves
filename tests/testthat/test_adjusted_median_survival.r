
set.seed(42)
sim_dat <- sim_confounded_surv(n=50, max_t=1.5)
sim_dat$group <- factor(sim_dat$group)

adj <- adjustedsurv(data=sim_dat,
                    variable="group",
                    ev_time="time",
                    event="event",
                    method="km",
                    conf_int=TRUE,
                    bootstrap=FALSE)

test_that("median surv, no boot", {
  adj_med <- adjusted_median_survival(adj)
  expect_equal(round(adj_med$median_surv, 4), c(0.4681, 0.6178))
})

test_that("median surv, with verbose", {
  expect_error(adjusted_median_survival(adj, verbose=FALSE), NA)
})
