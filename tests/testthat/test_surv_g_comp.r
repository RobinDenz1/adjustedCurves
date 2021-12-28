library(survival)
library(riskRegression)
library(pec)

set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_surv_n_50.Rds",
                               package="adjustedCurves"))
sim_dat$group <- as.factor(sim_dat$group)
sim_dat$event <- 1

mod <- glm(time ~ x1 + x2 + x4, data=sim_dat, family="gaussian")
class(mod) <- "custom"

# simplified version of predictProb.glm from the "pec" package,
# used only for testing here
custom_pred_fun <- function(object, newdata, times, ...) {
  N <- NROW(newdata)
  NT <- length(times)
  betax <- predict.glm(object, newdata=newdata, se.fit=FALSE)

  pred.matrix <- matrix(rep(times, N), byrow=TRUE, ncol=NT, nrow=N)
  p <- 1 - pnorm(pred.matrix - betax, mean=0, sd=sqrt(var(object$y)))
  return(p)
}

test_that("adjustedsurv, using predict_fun", {
  # using customly supplied fun
  adj <- suppressWarnings(adjustedsurv(data=sim_dat,
                                       variable="group",
                                       ev_time="time",
                                       event="event",
                                       method="direct",
                                       predict_fun=custom_pred_fun,
                                       outcome_model=mod,
                                       clean_data=FALSE))
  expect_s3_class(adj, "adjustedsurv")
  expect_true(is.numeric(adj$adjsurv$surv))
  expect_equal(levels(adj$adjsurv$group), levels(sim_dat$group))

  # using standard code
  class(mod) <- c("glm", "lm")
  adj_glm <- adjustedsurv(data=sim_dat,
                          variable="group",
                          ev_time="time",
                          event="event",
                          method="direct",
                          outcome_model=mod)
  expect_s3_class(adj_glm, "adjustedsurv")
  expect_true(is.numeric(adj_glm$adjsurv$surv))
  expect_equal(levels(adj_glm$adjsurv$group), levels(sim_dat$group))

  # should be equal, apart from the data and call parts
  adj_glm$data <- NULL
  adj_glm$call <- NULL
  adj$data <- NULL
  adj$call <- NULL
  expect_equal(adj, adj_glm)
})

test_that("adjustedsurv, using faulty S3 predict fun", {
  class(mod) <- "custom"

  expect_error(adjustedsurv(data=sim_dat,
                            variable="group",
                            ev_time="time",
                            event="event",
                            method="direct",
                            outcome_model=mod,
                            clean_data=FALSE),
               paste0("The following error occured using the default S3 ",
                      "predict method: 'Error in UseMethod(\"predict\"): ",
                      "no applicable method for 'predict' applied to an ",
                      "object of class \"custom\"\n' Specify a valid ",
                      "'predict_fun' or use a different model. See details."),
               fixed=TRUE)
})
