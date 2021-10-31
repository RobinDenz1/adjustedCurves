
adjsurv <- list(adjsurv=data.frame(time=seq(1, 100)))
class(adjsurv) <- "adjustedsurv"

not_adjsurv <- list()
class(not_adjsurv) <- "adjustedcif"

test_that("not an adjustedsurv object", {
  expect_error(adjustedCurves:::check_inputs_adj_rmst(adjsurv=not_adjsurv,
                                                      from=0,
                                                      to=1,
                                                      use_boot=FALSE),
               NULL)
})

test_that("from smaller 0", {
  expect_error(adjustedCurves:::check_inputs_adj_rmst(adjsurv=adjsurv,
                                                      from=-10,
                                                      to=1,
                                                      use_boot=FALSE),
               NULL)
})

test_that("from wrong format", {
  expect_error(adjustedCurves:::check_inputs_adj_rmst(adjsurv=adjsurv,
                                                      from="0",
                                                      to=1,
                                                      use_boot=FALSE),
               NULL)
})

test_that("from not smaller than to", {
  expect_error(adjustedCurves:::check_inputs_adj_rmst(adjsurv=adjsurv,
                                                      from=0,
                                                      to=0,
                                                      use_boot=FALSE),
               NULL)
})

test_that("no bootstrapping performed", {
  expect_warning(adjustedCurves:::check_inputs_adj_rmst(adjsurv=adjsurv,
                                                      from=0,
                                                      to=1,
                                                      use_boot=TRUE),
               NULL)
})
