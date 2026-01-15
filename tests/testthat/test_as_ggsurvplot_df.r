
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

test_that("works with multiple imputation", {

  set.seed(42)

  sim_dat <- readRDS(system.file("testdata",
                                 "d_sim_surv_n_50.Rds",
                                 package="adjustedCurves"))
  sim_dat$group <- as.factor(sim_dat$group)
  sim_dat$x1 <- ifelse(runif(n=nrow(sim_dat)) <= 0.7, sim_dat$x1, NA)

  # impute dataset
  imp <- suppressWarnings(mice::mice(sim_dat, m=3, method="pmm",
                                     printFlag=FALSE))

  # treatment model
  treat_mod <- with(imp, glm(group ~ x1 + x2 + x3 + x4 + x5 + x6,
                             family="binomial"))

  adj <- adjustedsurv(data=imp,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_km",
                      conf_int=FALSE,
                      treatment_model=treat_mod)

  out <- as_ggsurvplot_df(adj)

  expect_equal(colnames(out), c("time", "strata", "surv", "n.risk", "n.event"))
})
