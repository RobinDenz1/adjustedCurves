library(survival)

set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_surv_n_50.Rds",
                               package="adjustedCurves"))
sim_dat$group <- as.factor(sim_dat$group)
sim_dat$x2 <- ifelse(runif(n=nrow(sim_dat)) <= 0.5, sim_dat$x2, NA)

imp <- suppressWarnings(mice::mice(sim_dat, m=3, method="pmm",
                                   printFlag=FALSE))
treat_mod <- with(imp, glm(group ~ x2 + x3 + x5 + x6,
                           family="binomial"))

adj <- adjustedsurv(data=imp,
                    variable="group",
                    ev_time="time",
                    event="event",
                    method="iptw_km",
                    treatment_model=treat_mod,
                    bootstrap=TRUE,
                    n_boot=5)

# Rubins Rule applied to the results of the single imputed datasets
rubin_se <- function(est, se) {
  m <- length(est)
  sqrt(mean(se^2) + (1 + 1/m) * stats::var(est))
}

test_that("adjusted_rmst pools the bootstrap se using Rubins Rule", {
  pooled <- adjusted_rmst(adj, to=0.8, conf_int=TRUE)
  single <- lapply(adj$mids_analyses, adjusted_rmst, to=0.8, conf_int=TRUE)
  rmst <- vapply(single, function(x) x$rmst, numeric(2))
  se <- vapply(single, function(x) x$se, numeric(2))

  expect_true(all(apply(rmst, 1, stats::var) > 0))
  expect_equal(pooled$rmst, rowMeans(rmst))
  expect_equal(pooled$se, c(rubin_se(rmst[1, ], se[1, ]),
                            rubin_se(rmst[2, ], se[2, ])))
})

test_that("adjusted_rmtl pools the bootstrap se using Rubins Rule", {
  pooled <- adjusted_rmtl(adj, to=0.8, conf_int=TRUE)
  single <- lapply(adj$mids_analyses, adjusted_rmtl, to=0.8, conf_int=TRUE)
  rmtl <- vapply(single, function(x) x$rmtl, numeric(2))
  se <- vapply(single, function(x) x$se, numeric(2))

  expect_equal(pooled$se, c(rubin_se(rmtl[1, ], se[1, ]),
                            rubin_se(rmtl[2, ], se[2, ])))
})
