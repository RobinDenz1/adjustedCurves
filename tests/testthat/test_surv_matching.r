
set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_surv_n_100.Rds",
                               package="adjustedCurves"))
sim_dat$group <- factor(sim_dat$group)

# estimate propensity score
mod <- glm(group ~ x1 + x2 + x3 + x4 + x5 + x6, data=sim_dat,
           family="binomial")
ps_score <- mod$fitted.values

## Just check if function throws any errors
test_that("2 treatments, no conf_int, no boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="matching",
                      conf_int=FALSE,
                      treatment_model=mod)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("2 treatments, no conf_int, no boot, with ps_score", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="matching",
                      conf_int=FALSE,
                      bootstrap=FALSE,
                      n_boot=2,
                      treatment_model=ps_score)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})
