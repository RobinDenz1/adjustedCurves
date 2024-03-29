
set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_crisk_n_100.Rds",
                               package="adjustedCurves"))
sim_dat$group <- ifelse(sim_dat$group==0, "Control", "Treatment")
sim_dat$group <- as.factor(sim_dat$group)

mod <- glm(group ~ x1 + x2 + x5 + x6, data=sim_dat,
           family="binomial")

## Just check if function throws any errors
test_that("2 treatments, no conf_int, no boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="iptw_pseudo",
                     conf_int=FALSE,
                     treatment_model=mod,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("2 treatments, with conf_int, no boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="iptw_pseudo",
                     conf_int=TRUE,
                     treatment_model=mod,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("2 treatments, no conf_int, with boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="iptw_pseudo",
                     conf_int=FALSE,
                     bootstrap=TRUE,
                     n_boot=2,
                     treatment_model=mod,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("2 treatments, with conf_int, with boot", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="iptw_pseudo",
                     conf_int=TRUE,
                     bootstrap=TRUE,
                     n_boot=2,
                     treatment_model=mod,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("2 treatments, no conf_int, with WeightIt", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="iptw_pseudo",
                     conf_int=FALSE,
                     bootstrap=FALSE,
                     treatment_model=group ~ x1 + x2,
                     weight_method="ps",
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("2 treatments, no conf_int, with user-weights", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="iptw_pseudo",
                     conf_int=FALSE,
                     bootstrap=FALSE,
                     treatment_model=runif(n=50, min=1, max=2),
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("3 ways of iptw calculation are equal", {

  # calculate iptw curves in all 3 ways
  ps <- predict(mod, newdata=sim_dat, type="response")
  weights <- ifelse(sim_dat$group=="Treatment", 1/ps, 1/ (1-ps))

  adj_w <- adjustedcif(data=sim_dat,
                       variable="group",
                       ev_time="time",
                       event="event",
                       method="iptw_pseudo",
                       conf_int=TRUE,
                       bootstrap=FALSE,
                       treatment_model=weights,
                       stabilize=FALSE,
                       cause=1)$adj
  adj_glm <- adjustedcif(data=sim_dat,
                         variable="group",
                         ev_time="time",
                         event="event",
                         method="iptw_pseudo",
                         conf_int=TRUE,
                         bootstrap=FALSE,
                         treatment_model=mod,
                         stabilize=FALSE,
                         cause=1)$adj
  adj_weightit <- adjustedcif(data=sim_dat,
                              variable="group",
                              ev_time="time",
                              event="event",
                              method="iptw_pseudo",
                              conf_int=TRUE,
                              bootstrap=FALSE,
                              treatment_model=group ~ x1 + x2 + x5 + x6,
                              weight_method="ps",
                              stabilize=FALSE,
                              cause=1)$adj

  colnames(adj_w) <- paste0("w_", colnames(adj_w))
  colnames(adj_glm) <- paste0("glm_", colnames(adj_glm))
  colnames(adj_weightit) <- paste0("weightit_", colnames(adj_weightit))

  out <- cbind(adj_w, adj_glm, adj_weightit)
  out <- out[!is.na(out$glm_se),]

  tol <- 0.0001
  # all times equal
  expect_true(all.equal(out$w_time, out$glm_time, tolerance=tol))
  expect_true(all.equal(out$w_time, out$weightit_time, tolerance=tol))
  expect_true(all.equal(out$glm_time, out$weightit_time, tolerance=tol))

  # all surv equal
  expect_true(all.equal(out$w_cif, out$glm_cif, tolerance=tol))
  expect_true(all.equal(out$w_cif, out$weightit_cif, tolerance=tol))
  expect_true(all.equal(out$glm_cif, out$weightit_cif, tolerance=tol))

  # all groups
  expect_true(all.equal(out$w_group, out$glm_group, tolerance=tol))
  expect_true(all.equal(out$w_group, out$weightit_group, tolerance=tol))
  expect_true(all.equal(out$glm_group, out$weightit_group, tolerance=tol))

  # all se equal
  expect_true(all.equal(out$w_se, out$glm_se, tolerance=tol))
  # NOTE: On the devel version of R this throws an error saying there are
  #       NA values in out$glm_se, which is impossible since I remove all of
  #       those in line 148. I cannot replicate this bug on a local machine
  #       so I have no idea what causes this. Commented out until I can
  #       replicate & fix it
  #expect_true(all.equal(out$w_se, out$weightit_se, tolerance=tol))
  #expect_true(all.equal(out$glm_se, out$weightit_se, tolerance=tol))

})

set.seed(44522)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_crisk_n_150.Rds",
                               package="adjustedCurves"))
sim_dat$group[sim_dat$group==1] <- sample(c(1, 2),
                                        size=nrow(sim_dat[sim_dat$group==1, ]),
                                        replace=TRUE)
sim_dat$group <- as.factor(sim_dat$group)

mod <- quiet(nnet::multinom(group ~ x1 + x2 + x3 + x4 + x5 + x6, data=sim_dat))

test_that("> 2 treatments, no conf_int, no boot, no ...", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="iptw_pseudo",
                     conf_int=FALSE,
                     treatment_model=mod,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("> 2 treatments, with conf_int, no boot, no ...", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="iptw_pseudo",
                     conf_int=TRUE,
                     treatment_model=mod,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("> 2 treatments, no conf_int, with boot, no ...", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="iptw_pseudo",
                     conf_int=FALSE,
                     bootstrap=TRUE,
                     n_boot=2,
                     treatment_model=mod,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})

test_that("> 2 treatments, with conf_int, with boot, no ...", {
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="iptw_pseudo",
                     conf_int=TRUE,
                     bootstrap=TRUE,
                     n_boot=2,
                     treatment_model=mod,
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
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
  adj <- adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="iptw_pseudo",
                     conf_int=FALSE,
                     bootstrap=FALSE,
                     treatment_model=runif(n=100, min=1, max=2),
                     cause=1)
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))
})
