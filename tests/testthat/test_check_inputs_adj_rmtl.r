
adjsurv <- list(adjsurv=data.frame(time=seq(1, 100)))
class(adjsurv) <- "adjustedsurv"

not_adjsurv <- list()
class(not_adjsurv) <- "adjustedbla"

test_that("not an adjustedsurv/adjustedcif object", {
  expect_error(adjustedCurves:::check_inputs_adj_rmtl(adj=not_adjsurv,
                                                      from=0,
                                                      to=1,
                                                      use_boot=FALSE),
               NULL)
})

test_that("from smaller 0", {
  expect_error(adjustedCurves:::check_inputs_adj_rmtl(adj=adjsurv,
                                                      from=-10,
                                                      to=1,
                                                      use_boot=FALSE),
               "'from' and 'to' must be >= 0.")
})

test_that("from wrong format", {
  expect_error(adjustedCurves:::check_inputs_adj_rmtl(adj=adjsurv,
                                                      from="0",
                                                      to=1,
                                                      use_boot=FALSE),
               NULL)
})

test_that("from not smaller than to", {
  expect_error(adjustedCurves:::check_inputs_adj_rmtl(adj=adjsurv,
                                                      from=0,
                                                      to=0,
                                                      use_boot=FALSE),
               "'from' must be smaller than 'to'.")
})

test_that("no bootstrapping performed", {
  expect_warning(adjustedCurves:::check_inputs_adj_rmtl(adj=adjsurv,
                                                      from=0,
                                                      to=1,
                                                      use_boot=TRUE),
               paste0("Cannot use bootstrapped estimates because they ",
                      "were not estimated. Need 'bootstrap=TRUE' in ",
                      "'adjustedsurv'/'adjustedcif' function call."))
})

test_that("use_boot wrong format", {
  expect_error(adjustedCurves:::check_inputs_adj_rmst(adjsurv=adjsurv,
                                                      from=0,
                                                      to=1,
                                                      use_boot=1),
               "'use_boot' must be either TRUE or FALSE.")
})

test_that("no extrapolation allowed", {
  expect_error(adjustedCurves:::check_inputs_adj_rmst(adjsurv=adjsurv,
                                                      from=0,
                                                      to=200,
                                                      use_boot=FALSE),
               "'to' cannot be greater than the latest observed time.")
})

adjsurv <- list(adjsurv=data.frame(time=seq(1, 3)))
class(adjsurv) <- "adjustedsurv"

test_that("too little points in time", {
  expect_warning(adjustedCurves:::check_inputs_adj_rmst(adjsurv=adjsurv,
                                                        from=0,
                                                        to=1,
                                                        use_boot=FALSE),
                 paste0("Using only a few points in time might lead to ",
                        "biased estimates. Consider using a finer times ",
                        "grid in 'adjustedsurv'."))
})
