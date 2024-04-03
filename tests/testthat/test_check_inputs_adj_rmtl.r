
suppressMessages(requireNamespace("survival"))

adjsurv <- list(adjsurv=data.frame(time=seq(1, 100)))
class(adjsurv) <- "adjustedsurv"

not_adjsurv <- list()
class(not_adjsurv) <- "adjustedbla"

test_that("not an adjustedsurv/adjustedcif object", {
  expect_error(check_inputs_adj_rmtl(adj=not_adjsurv,
                                     from=0,
                                     to=1,
                                     conf_int=FALSE,
                                     contrast="none"),
               paste0("'adj' must be an 'adjustedsurv' object created ",
                      "using the 'adjustedsurv()' function or an ",
                      "'adjustedcif' object created using the ",
                      "'adjustedcif()' function."), fixed=TRUE)
})

test_that("from smaller 0", {
  expect_error(check_inputs_adj_rmtl(adj=adjsurv,
                                     from=-10,
                                     to=1,
                                     conf_int=FALSE,
                                     contrast="none"),
               "'from' and 'to' must be >= 0.")
})

test_that("from wrong format", {
  expect_error(check_inputs_adj_rmtl(adj=adjsurv,
                                     from="0",
                                     to=1,
                                     conf_int=FALSE,
                                     contrast="none"),
               "'from' and 'to' must be numbers (one for each argument).",
               fixed=TRUE)
})

test_that("from not smaller than to", {
  expect_error(check_inputs_adj_rmtl(adj=adjsurv,
                                     from=0,
                                     to=0,
                                     conf_int=FALSE,
                                     contrast="none"),
               "'from' must be smaller than 'to'.")
})

test_that("no bootstrapping performed", {
  expect_warning(check_inputs_adj_rmtl(adj=adjsurv,
                                       from=0,
                                       to=1,
                                       conf_int=TRUE,
                                       contrast="none"),
               paste0("Cannot use bootstrapped estimates because they ",
                      "were not estimated. Need 'bootstrap=TRUE' in ",
                      "'adjustedsurv'/'adjustedcif' function call."))
})

test_that("conf_int wrong format", {
  expect_error(check_inputs_adj_rmtl(adj=adjsurv,
                                     from=0,
                                     to=1,
                                     conf_int=1,
                                     contrast="none"),
               "'conf_int' must be either TRUE or FALSE.")
})

test_that("no extrapolation allowed", {
  expect_error(check_inputs_adj_rmtl(adj=adjsurv,
                                     from=0,
                                     to=200,
                                     conf_int=FALSE,
                                     contrast="none"),
               "'to' cannot be greater than the latest observed time.")
})

adjsurv <- list(adjsurv=data.frame(time=seq(1, 3)))
class(adjsurv) <- "adjustedsurv"

test_that("too little points in time", {
  expect_warning(check_inputs_adj_rmtl(adj=adjsurv,
                                       from=0,
                                       to=1,
                                       conf_int=FALSE,
                                       contrast="none"),
                 paste0("Using only a few points in time might lead to ",
                        "biased estimates. Consider using a finer times ",
                        "grid in 'adjustedsurv'."))
})
