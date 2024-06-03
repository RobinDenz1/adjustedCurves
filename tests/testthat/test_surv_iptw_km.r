
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
                      method="iptw_km",
                      conf_int=FALSE,
                      treatment_model=mod)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("2 treatments, with conf_int, no boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_km",
                      conf_int=TRUE,
                      treatment_model=mod)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("2 treatments, no conf_int, with boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_km",
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=2,
                      treatment_model=mod)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("2 treatments, with conf_int, with boot", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_km",
                      conf_int=TRUE,
                      bootstrap=TRUE,
                      n_boot=2,
                      treatment_model=mod)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("2 treatments, no conf_int, with WeightIt", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_km",
                      conf_int=FALSE,
                      bootstrap=FALSE,
                      treatment_model=group ~ x1 + x2,
                      weight_method="ps")
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("2 treatments, no conf_int, with WeightIt + additional args", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_km",
                      conf_int=FALSE,
                      bootstrap=FALSE,
                      treatment_model=group ~ x1 + x2,
                      weight_method="ps",
                      estimand="ATT")

  adj2 <- adjustedsurv(data=sim_dat,
                       variable="group",
                       ev_time="time",
                       event="event",
                       method="iptw_km",
                       conf_int=FALSE,
                       bootstrap=FALSE,
                       treatment_model=group ~ x1 + x2,
                       weight_method="ps",
                       estimand="ATE")

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
                      method="iptw_km",
                      conf_int=FALSE,
                      bootstrap=FALSE,
                      treatment_model=runif(n=50, min=1, max=2))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("3 ways of iptw calculation are equal", {

  # calculate iptw curves in all 3 ways
  ps <- predict(mod, newdata=sim_dat, type="response")
  weights <- ifelse(sim_dat$group=="Treatment", 1/ps, 1/ (1-ps))

  adj_w <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="iptw_km",
                        conf_int=TRUE,
                        bootstrap=FALSE,
                        treatment_model=weights,
                        stabilize=FALSE)$adj
  adj_glm <- adjustedsurv(data=sim_dat,
                          variable="group",
                          ev_time="time",
                          event="event",
                          method="iptw_km",
                          conf_int=TRUE,
                          bootstrap=FALSE,
                          treatment_model=mod,
                          stabilize=FALSE)$adj
  adj_weightit <- adjustedsurv(data=sim_dat,
                               variable="group",
                               ev_time="time",
                               event="event",
                               method="iptw_km",
                               conf_int=TRUE,
                               bootstrap=FALSE,
                               treatment_model=group ~ x1 + x2 + x3 +
                                 x4 + x5 + x6,
                               weight_method="ps",
                               stabilize=FALSE)$adj

  colnames(adj_w) <- paste0("w_", colnames(adj_w))
  colnames(adj_glm) <- paste0("glm_", colnames(adj_glm))
  colnames(adj_weightit) <- paste0("weightit_", colnames(adj_weightit))

  out <- cbind(adj_w, adj_glm, adj_weightit)

  tol <- 0.0001
  # all times equal
  expect_true(all.equal(out$w_time, out$glm_time, tolerance=tol))
  expect_true(all.equal(out$w_time, out$weightit_time, tolerance=tol))
  expect_true(all.equal(out$glm_time, out$weightit_time, tolerance=tol))

  # all surv equal
  expect_true(all.equal(out$w_surv, out$glm_surv, tolerance=tol))
  expect_true(all.equal(out$w_surv, out$weightit_surv, tolerance=tol))
  expect_true(all.equal(out$glm_surv, out$weightit_surv, tolerance=tol))

  # all groups
  expect_true(all.equal(out$w_group, out$glm_group, tolerance=tol))
  expect_true(all.equal(out$w_group, out$weightit_group, tolerance=tol))
  expect_true(all.equal(out$glm_group, out$weightit_group, tolerance=tol))

  # all se equal
  expect_true(all.equal(out$w_se, out$glm_se, tolerance=tol))
  expect_true(all.equal(out$w_se, out$weightit_se, tolerance=tol))
  expect_true(all.equal(out$glm_se, out$weightit_se, tolerance=tol))

})

test_that("use trim_quantiles", {
  adjtrim <- adjustedsurv(data=sim_dat,
                          variable="group",
                          ev_time="time",
                          event="event",
                          method="iptw_km",
                          conf_int=TRUE,
                          bootstrap=TRUE,
                          n_boot=2,
                          treatment_model=mod,
                          trim_quantiles=c(0.4, 0.6))
  expect_true(all(adjtrim$weights > 1.46 & adjtrim$weights < 1.74))
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
                      method="iptw_km",
                      conf_int=FALSE,
                      treatment_model=mod)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("> 2 treatments, with conf_int, no boot, no ...", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_km",
                      conf_int=TRUE,
                      treatment_model=mod)
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
                      treatment_model=mod)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("> 2 treatments, with conf_int, with boot, no ...", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_km",
                      conf_int=TRUE,
                      bootstrap=TRUE,
                      n_boot=2,
                      treatment_model=mod)
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

# DONT RUN: would require package dependency on "mlogit"
#test_that("> 2 treatments, no conf_int, with WeightIt", {
#  expect_error(adjustedCurves::adjustedsurv(data=sim_dat,
#                                            variable="group",
#                                            ev_time="time",
#                                            event="event",
#                                            method="iptw_km",
#                                            conf_int=FALSE,
#                                            bootstrap=FALSE,
#                                            treatment_model=group ~ x1 + x2,
#                                            weight_method="ps"), NA)
#})

test_that("> 2 treatments, no conf_int, with user-weights", {
  adj <- adjustedsurv(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw_km",
                      conf_int=FALSE,
                      bootstrap=FALSE,
                      treatment_model=runif(n=50, min=1, max=2))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adj$surv))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})
