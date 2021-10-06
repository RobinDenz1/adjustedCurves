
set.seed(42)

adjsurv <- list(boot_data=data.frame(a=1),
                adjsurv=data.frame(time=runif(20)))
class(adjsurv) <- "adjustedsurv"

adjsurv_no_boot <- list()
class(adjsurv_no_boot) <- "adjustedsurv"

not_adjsurv <- list()
class(not_adjsurv) <- "something"

test_that("not an adjustedsurv/adjustedcif object", {
  expect_error(adjustedCurves:::check_inputs_adj_test(adjsurv=not_adjsurv,
                                                      from=0,
                                                      to=1),
               NULL)
})

test_that("from smaller 0", {
  expect_error(adjustedCurves:::check_inputs_adj_test(adjsurv=adjsurv,
                                                      from=-10,
                                                      to=1),
               NULL)
})

test_that("from wrong format", {
  expect_error(adjustedCurves:::check_inputs_adj_test(adjsurv=adjsurv,
                                                      from="0",
                                                      to=1),
               NULL)
})

test_that("from not smaller than to", {
  expect_error(adjustedCurves:::check_inputs_adj_test(adjsurv=adjsurv,
                                                      from=0,
                                                      to=0),
               NULL)
})

test_that("no bootstrapping performed", {
  expect_error(adjustedCurves:::check_inputs_adj_test(adjsurv=adjsurv_no_boot,
                                                        from=0,
                                                        to=1),
                 NULL)
})

test_that("to larger than last time", {
  expect_error(adjustedCurves:::check_inputs_adj_test(adjsurv=adjsurv,
                                                      from=0,
                                                      to=2),
               NULL)
})

adjsurv <- list(boot_data=data.frame(a=1),
                adjsurv=data.frame(time=runif(5)))
class(adjsurv) <- "adjustedsurv"

test_that("not enough points in time", {
  expect_warning(adjustedCurves:::check_inputs_adj_test(adjsurv=adjsurv,
                                                        from=0,
                                                        to=0.2),
               NULL)
})
