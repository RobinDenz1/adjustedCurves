
set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_surv_n_50.Rds",
                               package="adjustedCurves"))
sim_dat$group <- ifelse(sim_dat$group==0, "Control", "Treatment")
sim_dat$group <- as.factor(sim_dat$group)

mod <- glm(group ~ x1 + x2 + x3 + x4 + x5 + x6, data=sim_dat,
           family="binomial")

## Just check if function throws any errors
test_that("2 treatments, no conf_int, no boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_cens",
                      conf_int=FALSE,
                      treatment_model=mod,
                      censoring_model= ~ x2 + x4)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("2 treatments, no conf_int, with boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_cens",
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=2,
                      treatment_model=mod,
                      censoring_model= ~ x2 + x4)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("2 treatments, no conf_int, with WeightIt", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_cens",
                      conf_int=FALSE,
                      bootstrap=FALSE,
                      treatment_model=group ~ x1 + x2,
                      weight_method="ps",
                      censoring_model= ~ x2 + x4)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("2 treatments, no conf_int, with WeightIt + additional args", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_cens",
                      conf_int=FALSE,
                      bootstrap=FALSE,
                      treatment_model=group ~ x1 + x2,
                      weight_method="ps",
                      estimand="ATT",
                      censoring_model= ~ x2 + x4)

  adj2 <- adjustedsurv(data=sim_dat,
                       variable="group",
                       ev_time="time",
                       event="event",
                       method="iptw_cens",
                       conf_int=FALSE,
                       bootstrap=FALSE,
                       treatment_model=group ~ x1 + x2,
                       weight_method="ps",
                       estimand="ATE",
                       censoring_model= ~ x2 + x4)

  adj$call <- NULL
  adj2$call <- NULL

  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
  expect_true(!all(adj$adj$surv==adj2$adj$surv))
})

test_that("2 treatments, no conf_int, with user-weights", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_cens",
                      conf_int=FALSE,
                      bootstrap=FALSE,
                      treatment_model=runif(n=50, min=1, max=2),
                      censoring_model= ~ x2 + x4)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("use trim_quantiles", {
  adjtrim <- adjustedsurv(data=sim_dat,
                          variable="group",
                          ev_time="time",
                          event="event",
                          method="iptw_cens",
                          conf_int=FALSE,
                          bootstrap=TRUE,
                          n_boot=2,
                          treatment_model=mod,
                          trim_quantiles=c(0.4, 0.6),
                          censoring_model= ~ x2 + x4)
  expect_true(all(adjtrim$weights > 1.46 & adjtrim$weights < 1.74))
})

test_that("changing the Cox model specification", {
  adj2 <- adjustedsurv(data=sim_dat,
                       variable="group",
                       ev_time="time",
                       event="event",
                       method="iptw_cens",
                       conf_int=FALSE,
                       bootstrap=TRUE,
                       n_boot=2,
                       treatment_model=mod,
                       censoring_model= ~ x4,
                       ties="breslow",
                       coxph_control=list(outer.max=100))

  expect_equal(adj2$censoring_models[[1]]$method, "breslow")
})

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_surv_n_50.Rds",
                               package="adjustedCurves"))
sim_dat$group[sim_dat$group==1] <- sample(c(1, 2),
                                        size=nrow(sim_dat[sim_dat$group==1, ]),
                                        replace=TRUE)
sim_dat$group <- as.factor(sim_dat$group)

mod <- quiet(nnet::multinom(group ~ x1 + x2 + x3 + x4 + x5 + x6, data=sim_dat))

test_that("> 2 treatments, no conf_int, no boot, no ...", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_cens",
                      conf_int=FALSE,
                      treatment_model=mod,
                      censoring_model= ~ x4)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("> 2 treatments, no conf_int, with boot, no ...", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_km",
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=2,
                      treatment_model=mod,
                      censoring_model= ~ x4)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("> 2 treatments, no conf_int, with user-weights", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_cens",
                      conf_int=FALSE,
                      bootstrap=FALSE,
                      treatment_model=runif(n=50, min=1, max=2),
                      censoring_model= ~ x4)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})
