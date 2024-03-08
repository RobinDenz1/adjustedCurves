
set.seed(42)

## 2 treatments
sim_dat <- readRDS(system.file("testdata",
                               "d_sim_surv_n_20.Rds",
                               package="adjustedCurves"))
sim_dat$group <- as.factor(sim_dat$group)

sim_dat_tibble <- dplyr::tibble(sim_dat)
sim_dat$group_num <- sim_dat$group
sim_dat$group <- as.factor(sim_dat$group)

# treatment model
treat_mod <- glm(group ~ x1 + x2 + x3 + x4 + x5 + x6, data=sim_dat,
                 family="binomial")

## Just check if function throws any errors
test_that("2 treatments, using glm", {
  weights <- get_iptw_weights(data=sim_dat,
                              treatment_model=treat_mod,
                              weight_method="ps",
                              variable="group",
                              stabilize=TRUE,
                              trim=FALSE,
                              trim_q=FALSE)
  expect_true(is.numeric(weights))
  expect_true(all(weights > 0))
  expect_true(length(weights)==nrow(sim_dat))
})

test_that("2 treatments, using weightit", {
  weights <- get_iptw_weights(data=sim_dat,
                              treatment_model=group ~ x2 + x3,
                              weight_method="ps",
                              variable="group",
                              stabilize=TRUE,
                              trim=FALSE,
                              trim_q=FALSE)
  expect_true(is.numeric(weights))
  expect_true(all(weights > 0))
  expect_true(length(weights)==nrow(sim_dat))
})

test_that("2 treatments, using glm + trim", {
  weights <- get_iptw_weights(data=sim_dat,
                              treatment_model=treat_mod,
                              weight_method="ps",
                              variable="group",
                              stabilize=TRUE,
                              trim=3,
                              trim_q=FALSE)
  expect_true(is.numeric(weights))
  expect_true(all(weights > 0))
  expect_true(length(weights)==nrow(sim_dat))
})

test_that("2 treatments, using weightit + trim", {
  weights <- get_iptw_weights(data=sim_dat,
                              treatment_model=group ~ x2 + x3,
                              weight_method="ps",
                              variable="group",
                              stabilize=TRUE,
                              trim=3,
                              trim_q=FALSE)
  expect_true(is.numeric(weights))
  expect_true(all(weights > 0))
  expect_true(length(weights)==nrow(sim_dat))
})

test_that("2 treatments, not using stabilize", {
  weights <- get_iptw_weights(data=sim_dat,
                              treatment_model=treat_mod,
                              weight_method="ps",
                              variable="group",
                              stabilize=TRUE,
                              trim=3,
                              trim_q=FALSE)
  expect_true(is.numeric(weights))
  expect_true(all(weights > 0))
  expect_true(length(weights)==nrow(sim_dat))
})

mira_mod <- list(A="stuff")
class(mira_mod) <- "mira"

test_that("error with mira", {
  expect_error(get_iptw_weights(data=sim_dat,
                                treatment_model=mira_mod,
                                weight_method="ps",
                                variable="group",
                                stabilize=TRUE,
                                trim=3,
                                trim_q=FALSE),
                paste0("If a 'mira' object is used in the 'treatment_model' ",
                       "argument the 'data' argument must be a 'mids' ",
                       "object, not a data.frame."))
})

class(mira_mod) <- "unsupported"

test_that("error unsupported input", {
  expect_error(get_iptw_weights(data=sim_dat,
                                treatment_model=mira_mod,
                                weight_method="ps",
                                variable="group",
                                stabilize=TRUE,
                                trim=3,
                                trim_q=FALSE),
               "Unsuported input: 'unsupported'. See documentation.")
})

## 3 treatments
sim_dat$group2 <- 0
sim_dat$group2[sim_dat$group_num==1] <-
  sample(c(1, 2), size=nrow(sim_dat[sim_dat$group_num==1, ]), replace=TRUE)
sim_dat$group2 <- ifelse(sim_dat$group2==1, "Placebo",
                         ifelse(sim_dat$group2==2, "Chemo", "OP"))
sim_dat$group2 <- factor(sim_dat$group2)

treat_mod <- quiet(nnet::multinom(group2 ~ x1 + x2 + x4, data=sim_dat))

## Just check if function throws any errors
test_that("3 treatments, using multinom", {
  weights <- get_iptw_weights(data=sim_dat,
                              treatment_model=treat_mod,
                              weight_method="ps",
                              variable="group2",
                              stabilize=TRUE,
                              trim=FALSE,
                              trim_q=FALSE)
  expect_true(is.numeric(weights))
  expect_true(all(weights > 0))
  expect_true(length(weights)==nrow(sim_dat))
})

test_that("3 treatments, using weightit", {
  weights <- get_iptw_weights(data=sim_dat,
                              treatment_model=group ~ x2 + x3,
                              weight_method="ps",
                              variable="group2",
                              stabilize=TRUE,
                              trim=FALSE,
                              trim_q=FALSE)
  expect_true(is.numeric(weights))
  expect_true(all(weights > 0))
  expect_true(length(weights)==nrow(sim_dat))
})

test_that("3 treatments, using multinom + trim", {
  weights <- get_iptw_weights(data=sim_dat,
                              treatment_model=treat_mod,
                              weight_method="ps",
                              variable="group2",
                              stabilize=TRUE,
                              trim=3,
                              trim_q=FALSE)
  expect_true(is.numeric(weights))
  expect_true(all(weights > 0))
  expect_true(length(weights)==nrow(sim_dat))
})

test_that("3 treatments, using weightit + trim", {
  weights <- get_iptw_weights(data=sim_dat,
                              treatment_model=group ~ x2 + x3,
                              weight_method="ps",
                              variable="group2",
                              stabilize=TRUE,
                              trim=3,
                              trim_q=FALSE)
  expect_true(is.numeric(weights))
  expect_true(all(weights > 0))
  expect_true(length(weights)==nrow(sim_dat))
})

test_that("3 treatments, not using stabilize", {
  weights <- get_iptw_weights(data=sim_dat,
                              treatment_model=treat_mod,
                              weight_method="ps",
                              variable="group2",
                              stabilize=TRUE,
                              trim=3,
                              trim_q=FALSE)
  expect_true(is.numeric(weights))
  expect_true(all(weights > 0))
  expect_true(length(weights)==nrow(sim_dat))
})
