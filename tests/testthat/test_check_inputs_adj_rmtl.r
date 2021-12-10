
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
               NULL)
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
               NULL)
})

test_that("no bootstrapping performed", {
  expect_warning(adjustedCurves:::check_inputs_adj_rmtl(adj=adjsurv,
                                                      from=0,
                                                      to=1,
                                                      use_boot=TRUE),
               NULL)
})

test_that("use_boot wrong format", {
  expect_error(adjustedCurves:::check_inputs_adj_rmst(adjsurv=adjsurv,
                                                      from=0,
                                                      to=1,
                                                      use_boot=1),
               NULL)
})

test_that("no extrapolation allowed", {
  expect_error(adjustedCurves:::check_inputs_adj_rmst(adjsurv=adjsurv,
                                                      from=0,
                                                      to=200,
                                                      use_boot=FALSE),
               NULL)
})

adjsurv <- list(adjsurv=data.frame(time=seq(1, 3)))
class(adjsurv) <- "adjustedsurv"

test_that("too little points in time", {
  expect_warning(adjustedCurves:::check_inputs_adj_rmst(adjsurv=adjsurv,
                                                        from=0,
                                                        to=1,
                                                        use_boot=FALSE),
                 NULL)
})

