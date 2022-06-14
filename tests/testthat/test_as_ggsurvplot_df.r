
library(survival)
library(prodlim)

set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_surv_n_20.Rds",
                               package="adjustedCurves"))
sim_dat$group <- as.factor(sim_dat$group)

adj <- adjustedsurv(data=sim_dat,
                    variable="group",
                    event="event",
                    ev_time="time",
                    method="iptw_km",
                    treatment_model=group ~ x1)

test_that("as_ggsurvplot_df works as expected", {
  df <- as_ggsurvplot_df(adj)

  expect_true(is.data.frame(df))
  expect_equal(colnames(df), c("time", "surv", "strata", "n.risk", "n.event"))
})

test_that("as_ggsurvplot_df wrong input", {
  expect_error(as_ggsurvplot_df("WRONG"),
               paste0("This function can only be used with 'adjustedsurv'",
                      " objects, created using the 'adjustedsurv' function."))
})

