library(prodlim)

set.seed(41)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_crisk_n_100.Rds",
                               package="adjustedCurves"))
sim_dat$group <- as.factor(sim_dat$group)

sim_dat$x1 <- ifelse(stats::runif(n=100)<0.5, NA, sim_dat$x1)

mids <- mice::mice(data=sim_dat, m=3, method="pmm", printFlag=FALSE)

test_that("no additional args", {
  mod <- FGR_MI(mids=mids, formula=Hist(time, event) ~ x1 + x2)
  expect_s3_class(mod, "mira")
})

test_that("with additional args", {
  mod <- FGR_MI(mids=mids, formula=Hist(time, event) ~ x1 + x2, cause=1)
  expect_s3_class(mod, "mira")
})

test_that("trying to use data", {
  expect_error(FGR_MI(mids=mids,
                      formula=Hist(time, event) ~ x1 + x2,
                      cause=1,
                      data=sim_dat),
               paste0("The 'data' argument cannot be used here. ",
                      "It is replaced by the'mids' argument."))
})
