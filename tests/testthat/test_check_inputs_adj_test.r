
suppressMessages(requireNamespace("survival"))

set.seed(42)

adjsurv <- list(boot_data=data.frame(a=1),
                adj=data.frame(time=runif(20)),
                adjcif=data.frame(time=runif(20)))
class(adjsurv) <- "adjustedsurv"

adjsurv_no_boot <- list()
class(adjsurv_no_boot) <- "adjustedsurv"

not_adjsurv <- list()
class(not_adjsurv) <- "something"

test_that("not an adjustedsurv/adjustedcif object", {
  expect_error(check_inputs_adj_test(adj=not_adjsurv, from=0, to=1),
               paste0("'adj' must be an 'adjustedsurv' or ",
                      "'adjustedcif' object,created using the ",
                      "adjustedsurv or adjustedcif function."))
})

test_that("from smaller 0", {
  expect_error(check_inputs_adj_test(adj=adjsurv, from=-10, to=1),
               "'from' must be a number >= 0.")
})

test_that("from wrong format", {
  expect_error(check_inputs_adj_test(adj=adjsurv, from="0", to=1),
               "'from' must be a single number >= 0.")
})

test_that("to wrong format", {
  expect_error(check_inputs_adj_test(adj=adjsurv, from=0, to="1"),
               "'to' must be a single number >= 0.")
})

test_that("adjustedsurv, from not smaller than to", {
  expect_error(check_inputs_adj_test(adj=adjsurv, from=0, to=0),
               "'to' must be greater than 'from'.")
})

test_that("adjustedsurv, no bootstrapping performed", {
  expect_error(check_inputs_adj_test(adj=adjsurv_no_boot, from=0, to=1),
               paste0("Can only perform a significance test if bootstrapping ",
                      "was performed (bootstrap=TRUE in ",
                      "adjustedsurv/adjustedcif call)."),
               fixed=TRUE)
})

test_that("adjustedsurv, to larger than last time", {
  expect_error(check_inputs_adj_test(adj=adjsurv, from=0, to=2),
               "'to' cannot be greater than the latest observed time.")
})

class(adjsurv) <- "adjustedcif"

test_that("adjustedcif, to larger than last time", {
  expect_error(check_inputs_adj_test(adj=adjsurv, from=0, to=2),
               "'to' cannot be greater than the latest observed time.")
})

test_that("adjustedsurv, no bootstrapping performed", {
  expect_error(check_inputs_adj_test(adj=adjsurv_no_boot, from=0, to=1),
               paste0("Can only perform a significance test if bootstrapping ",
                      "was performed (bootstrap=TRUE in ",
                      "adjustedsurv/adjustedcif call)."),
               fixed=TRUE)
})

class(adjsurv_no_boot) <- "adjustedcif"

test_that("adjustedcif, no bootstrapping performed", {
  expect_error(check_inputs_adj_test(adj=adjsurv_no_boot, from=0, to=1),
               paste0("Can only perform a significance test if bootstrapping ",
                      "was performed (bootstrap=TRUE in ",
                      "adjustedsurv/adjustedcif call)."),
               fixed=TRUE)
})

adjsurv <- list(boot_data=data.frame(a=1),
                adj=data.frame(time=runif(5)),
                adjcif=data.frame(time=runif(5)))
class(adjsurv) <- "adjustedsurv"

test_that("adjustedsurv, not enough points in time", {
  expect_warning(check_inputs_adj_test(adj=adjsurv, from=0, to=0.2),
               paste0("Using only a few points in time might lead to ",
                      "biased estimates. Consider using a finer times ",
                      "grid in 'adjustedsurv' or 'adjustedcif'."))
})

class(adjsurv) <- "adjustedcif"

test_that("adjustedcif, not enough points in time", {
  expect_warning(check_inputs_adj_test(adj=adjsurv, from=0, to=0.2),
                 paste0("Using only a few points in time might lead to ",
                        "biased estimates. Consider using a finer times ",
                        "grid in 'adjustedsurv' or 'adjustedcif'."))
})
