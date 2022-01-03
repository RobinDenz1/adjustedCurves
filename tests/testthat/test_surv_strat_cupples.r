library(survival)

set.seed(42)

sim_dat <- sim_confounded_surv(n=200, list(x1=c("rbinom", 1, 0.5),
                                           x2=c("rbinom", 1, 0.5),
                                           x3=c("rbinom", 1, 0.5)),
                               treatment_betas=c(x1=0.2, x2=0.4, x3=0.4),
                               outcome_betas=c(x1=1.1, x2=0, x3=1.1),
                               group_beta=0)
sim_dat$group2 <- factor(sim_dat$group)
reference <- sim_dat[, c("x1", "x3")]

## Just check if function throws any errors
test_that("2 treatments, one confounder, no boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group2",
                      ev_time="time",
                      event="event",
                      method="strat_cupples",
                      adjust_vars="x1")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group2))
})

test_that("2 treatments, two confounders, no boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group2",
                      ev_time="time",
                      event="event",
                      method="strat_cupples",
                      adjust_vars=c("x1", "x3"))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group2))
})

test_that("2 treatments, two confounders, with boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group2",
                      ev_time="time",
                      event="event",
                      method="strat_cupples",
                      adjust_vars=c("x1", "x3"),
                      bootstrap=TRUE,
                      n_boot=2)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group2))
})

test_that("2 treatments, two confounders, no boot, using reference data", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group2",
                      ev_time="time",
                      event="event",
                      method="strat_cupples",
                      adjust_vars=c("x1", "x3"),
                      reference=reference)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group2))
})

## more than two treatments
sim_dat$group[sim_dat$group==1] <- sample(c(1, 2),
                                        size=nrow(sim_dat[sim_dat$group==1, ]),
                                        replace=TRUE)
sim_dat$group <- factor(sim_dat$group)

test_that("> 2 treatments, one confounder, no boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="strat_cupples",
                      adjust_vars="x1")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("> 2 treatments, two confounders, no boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="strat_cupples",
                      adjust_vars=c("x1", "x3"))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("> 2 treatments, two confounders, with boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="strat_cupples",
                      adjust_vars=c("x1", "x3"),
                      bootstrap=TRUE,
                      n_boot=2)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})

test_that("> 2 treatments, two confounders, no boot, using reference data", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="strat_cupples",
                      adjust_vars=c("x1", "x3"),
                      reference=reference)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))
})
