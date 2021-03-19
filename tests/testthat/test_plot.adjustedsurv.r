library(survival)

sim_dat <- sim_confounded_surv(n=200)
sim_dat$group <- factor(sim_dat$group)

adj <- adjustedsurv(data=sim_dat,
                    variable="group",
                    ev_time="time",
                    event="event",
                    method="km",
                    conf_int=T,
                    bootstrap=T,
                    n_boot=10)

test_that("plot, no conf_int", {
  expect_error(plot(adj), NA)
})

test_that("plot, with conf_int", {
  expect_error(plot(adj, draw_ci=T), NA)
})

test_that("plot, no conf_int, using boot", {
  expect_error(plot(adj, draw_ci=T, use_boot=T), NA)
})

test_that("plot, using iso_reg", {
  expect_error(plot(adj, iso_reg=T), NA)
})
