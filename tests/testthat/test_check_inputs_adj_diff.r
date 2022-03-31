
suppressMessages(requireNamespace("survival"))

set.seed(42)

sim_dat <- readRDS(system.file("testdata",
                               "d_sim_surv_n_100.Rds",
                               package="adjustedCurves"))
sim_dat$group <- factor(sim_dat$group)

adj <- adjustedsurv(data=sim_dat,
                    variable="group",
                    ev_time="time",
                    event="event",
                    method="km",
                    conf_int=TRUE,
                    bootstrap=TRUE,
                    n_boot=2)

adj2 <- adj
adj2$adjcif <- adj$adjsurv
class(adj2) <- "adjustedcif"

test_that("not an adjustedsurv/adjustedcif object", {
  expect_error(check_inputs_adj_diff(adj="a",
                                     group_1="0",
                                     group_2="1",
                                     conf_int=TRUE,
                                     use_boot=FALSE),
               paste0("'x' must be an 'adjustedsurv' object created using",
                      " the adjustedsurv function or an 'adjustedcif' object ",
                      "created using the adjustedcif function."),
               fixed=TRUE)
})

test_that("wrong group_1", {
  expect_error(check_inputs_adj_diff(adj=adj,
                                     group_1=c("0", "1"),
                                     group_2="1",
                                     conf_int=TRUE,
                                     use_boot=FALSE),
               paste0("'group_1' has to be a single character vector, ",
                      "specifying one of the treatment groups in 'variable'."),
               fixed=TRUE)
})

test_that("wrong group_2", {
  expect_error(check_inputs_adj_diff(adj=adj,
                                     group_1="0",
                                     group_2=c("0", "1"),
                                     conf_int=TRUE,
                                     use_boot=FALSE),
               paste0("'group_2' has to be a single character vector, ",
                      "specifying one of the treatment groups in 'variable'."),
               fixed=TRUE)
})

test_that("group_1 not in adjustedsurv", {
  expect_error(check_inputs_adj_diff(adj=adj,
                                     group_1="2",
                                     group_2="0",
                                     conf_int=TRUE,
                                     use_boot=FALSE),
               paste0("2 is not a valid group in 'variable'."),
               fixed=TRUE)
})

test_that("group_1 not in adjustedcif", {
  expect_error(check_inputs_adj_diff(adj=adj2,
                                     group_1="2",
                                     group_2="0",
                                     conf_int=TRUE,
                                     use_boot=FALSE),
               paste0("2 is not a valid group in 'variable'."),
               fixed=TRUE)
})

test_that("group_2 not in adjustedsurv", {
  expect_error(check_inputs_adj_diff(adj=adj,
                                     group_1="0",
                                     group_2="2",
                                     conf_int=TRUE,
                                     use_boot=FALSE),
               paste0("2 is not a valid group in 'variable'."),
               fixed=TRUE)
})

test_that("group_2 not in adjustedcif", {
  expect_error(check_inputs_adj_diff(adj=adj2,
                                     group_1="0",
                                     group_2="2",
                                     conf_int=TRUE,
                                     use_boot=FALSE),
               paste0("2 is not a valid group in 'variable'."),
               fixed=TRUE)
})

test_that("group_1 same as group_2", {
  expect_error(check_inputs_adj_diff(adj=adj,
                                     group_1="0",
                                     group_2="0",
                                     conf_int=TRUE,
                                     use_boot=FALSE),
               paste0("'group_1' and 'group_2' may not be equal."),
               fixed=TRUE)
})

adj$boot_adjsurv <- NULL

test_that("using use_boot without boot", {
  expect_error(check_inputs_adj_diff(adj=adj,
                                     group_1="0",
                                     group_2="1",
                                     conf_int=TRUE,
                                     use_boot=TRUE),
               paste0("Bootstrapped estimates can only be calculated if ",
                      "'bootstrap=TRUE' was used in the original ",
                      "adjustedsurv or adjustedcif function call."),
               fixed=TRUE)
})

adj$adjsurv <- dplyr::select(adj$adjsurv, c("time", "group", "surv"))

test_that("using conf_int without approximate stuff surv", {
  expect_error(check_inputs_adj_diff(adj=adj,
                                     group_1="0",
                                     group_2="1",
                                     conf_int=TRUE,
                                     use_boot=FALSE),
               paste0("There are no approximate standard error ",
                      "calculations to use. Either set 'use_boot=TRUE' or ",
                      "rerun the adjustedsurv function with ",
                      "'conf_int=TRUE' if possible."),
               fixed=TRUE)
})

adj2$adjcif <- dplyr::select(adj2$adjcif, c("time", "group", "surv"))

test_that("using conf_int without approximate stuff surv", {
  expect_error(check_inputs_adj_diff(adj=adj2,
                                     group_1="0",
                                     group_2="1",
                                     conf_int=TRUE,
                                     use_boot=FALSE),
               paste0("There are no approximate standard error ",
                      "calculations to use. Either set 'use_boot=TRUE' or ",
                      "rerun the adjustedcif function with ",
                      "'conf_int=TRUE' if possible."),
               fixed=TRUE)
})
