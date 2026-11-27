
set.seed(42)

sim_dat <- sim_confounded_surv(n=200, list(x1=c("rbinom", 1, 0.5),
                                           x2=c("rbinom", 1, 0.5),
                                           x3=c("rbinom", 1, 0.5)),
                               treatment_betas=c(x1=0.2, x2=0.4, x3=0.4),
                               outcome_betas=c(x1=1.1, x2=0, x3=1.1),
                               group_beta=0)
sim_dat$group2 <- factor(sim_dat$group)


## Just check if function throws any errors
test_that("2 treatments, one confounder", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group2",
                      ev_time="time",
                      event="event",
                      method="strat_nieto",
                      adjust_vars="x1")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group2))
})

test_that("2 treatments, one confounder, with conf_int", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group2",
                      ev_time="time",
                      event="event",
                      method="strat_nieto",
                      adjust_vars="x1",
                      conf_int=TRUE)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group2))
})

test_that("se is the standard error of the survival probability", {
  dat_one <- sim_dat
  dat_one$const <- factor("all")
  times <- c(0.2, 0.4, 0.6)

  adj_nieto <- adjustedsurv(data=dat_one,
                            variable="group2",
                            ev_time="time",
                            event="event",
                            method="strat_nieto",
                            adjust_vars="const",
                            conf_int=TRUE,
                            times=times)
  adj_km <- adjustedsurv(data=dat_one,
                         variable="group2",
                         ev_time="time",
                         event="event",
                         method="km",
                         conf_int=TRUE,
                         times=times)

  expect_equal(adj_nieto$adj$surv, adj_km$adj$surv)
  expect_equal(adj_nieto$adj$se, adj_km$adj$se)
  expect_equal(adj_nieto$adj$ci_lower, adj_km$adj$ci_lower)
  expect_equal(adj_nieto$adj$ci_upper, adj_km$adj$ci_upper)

  diff_nieto <- adjusted_curve_diff(adj_nieto, times=times, conf_int=TRUE)
  diff_km <- adjusted_curve_diff(adj_km, times=times, conf_int=TRUE)
  expect_equal(diff_nieto$se, diff_km$se)
  expect_equal(diff_nieto$p_value, diff_km$p_value)
})

test_that("2 treatments, one confounder, with boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group2",
                      ev_time="time",
                      event="event",
                      method="strat_nieto",
                      adjust_vars="x1",
                      bootstrap=TRUE,
                      n_boot=2)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group2))
})

test_that("2 treatments, two confounders", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group2",
                      ev_time="time",
                      event="event",
                      method="strat_nieto",
                      adjust_vars=c("x1", "x3"))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group2))
})

test_that("2 treatments, two confounders, with conf_int", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group2",
                      ev_time="time",
                      event="event",
                      method="strat_nieto",
                      adjust_vars=c("x1", "x3"),
                      conf_int=TRUE)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group2))
})

test_that("2 treatments, two confounders, with boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group2",
                      ev_time="time",
                      event="event",
                      method="strat_nieto",
                      adjust_vars=c("x1", "x3"),
                      bootstrap=TRUE,
                      n_boot=2)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group2))
})

## more than two treatments
sim_dat$group[sim_dat$group==1] <- sample(c(1, 2),
                                        size=nrow(sim_dat[sim_dat$group==1, ]),
                                        replace=TRUE)
sim_dat$group <- factor(sim_dat$group)

test_that("> 2 treatments, one confounder", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="strat_nieto",
                      adjust_vars="x1")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("> 2 treatments, one confounder, with conf_int", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="strat_nieto",
                      adjust_vars="x1",
                      conf_int=TRUE)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("> 2 treatments, one confounder, with boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="strat_nieto",
                      adjust_vars="x1",
                      bootstrap=TRUE,
                      n_boot=2)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("> 2 treatments, two confounders", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="strat_nieto",
                      adjust_vars=c("x1", "x3"))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("> 2 treatments, two confounders, with conf_int", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="strat_nieto",
                      adjust_vars=c("x1", "x3"),
                      conf_int=TRUE)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("> 2 treatments, two confounders, with boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="strat_nieto",
                      adjust_vars=c("x1", "x3"),
                      bootstrap=TRUE,
                      n_boot=2)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})
