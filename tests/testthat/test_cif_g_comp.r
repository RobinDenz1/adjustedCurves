library(survival)
library(riskRegression)
library(pec)
library(prodlim)

set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_crisk_n_50.Rds",
                               package="adjustedCurves"))
sim_dat$group <- as.factor(sim_dat$group)

mod <- FGR(Hist(time, event) ~ group + x1 + x2, data=sim_dat, cause=1)
class(mod) <- "custom"

custom_pred_fun <- getFromNamespace("predict.FGR", "riskRegression")

test_that("adjustedcif, using predict_fun", {
  # using customly supplied fun
  adj <- quiet(adjustedcif(data=sim_dat,
                     variable="group",
                     ev_time="time",
                     event="event",
                     method="direct",
                     predict_fun=custom_pred_fun,
                     outcome_model=mod,
                     clean_data=FALSE,
                     cause=1))
  expect_s3_class(adj, "adjustedcif")
  expect_true(is.numeric(adj$adj$cif))
  expect_equal(levels(adj$adj$group), levels(sim_dat$group))

  # using standard code
  class(mod) <- "FGR"
  adj_fgr <- adjustedcif(data=sim_dat,
                         variable="group",
                         ev_time="time",
                         event="event",
                         method="direct",
                         outcome_model=mod,
                         cause=1)
  expect_s3_class(adj_fgr, "adjustedcif")
  expect_true(is.numeric(adj_fgr$adj$cif))
  expect_equal(levels(adj_fgr$adj$group), levels(sim_dat$group))

  # should be equal, apart from the data and call parts
  adj_fgr$data <- NULL
  adj_fgr$call <- NULL
  adj$data <- NULL
  adj$call <- NULL
  expect_equal(adj, adj_fgr)
})

test_that("adjustedcif, using faulty S3 predict fun", {
  class(mod) <- "custom"

  expect_error(adjustedcif(data=sim_dat,
                           variable="group",
                           ev_time="time",
                           event="event",
                           method="direct",
                           outcome_model=mod,
                           clean_data=FALSE,
                           cause=1),
               paste0("The following error occured using the default S3 ",
                      "predict method: 'Error in UseMethod(\"predict\"): ",
                      "no applicable method for 'predict' applied to an ",
                      "object of class \"custom\"\n' Specify a valid ",
                      "'predict_fun' or use a different model. See details."),
               fixed=TRUE)
})
