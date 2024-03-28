
set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_crisk_n_50.Rds",
                               package="adjustedCurves"))
sim_dat$group <- as.factor(sim_dat$group)

# estimate propensity score
mod <- glm(group ~ x1 + x2 + x3 + x4 + x5 + x6, data=sim_dat,
           family="binomial")
ps_score <- mod$fitted.values

test_that("2 treatments, no conf_int, no boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="matching",
                     conf_int=FALSE,
                     treatment_model=mod,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
})

test_that("2 treatments, no conf_int, no boot, with ps_score", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="matching",
                     conf_int=FALSE,
                     bootstrap=FALSE,
                     n_boot=2,
                     treatment_model=ps_score,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
})
