
set.seed(42)

values <- stats::rnorm(100)
weights <- stats::runif(100)

test_that("miller", {
  results <- weighted.var.se(x=values, w=weights, se_method="miller")
  expect_true(!is.na(results))
  expect_length(results, 1)
  expect_gt(results, 0)
})

test_that("galloway", {
  results <- weighted.var.se(x=values, w=weights, se_method="galloway")
  expect_true(!is.na(results))
  expect_length(results, 1)
  expect_gt(results, 0)
})

test_that("cochrane", {
  results <- weighted.var.se(x=values, w=weights, se_method="cochrane")
  expect_true(!is.na(results))
  expect_length(results, 1)
  expect_gt(results, 0)
})

test_that("simple", {
  results <- weighted.var.se(x=values, w=weights, se_method="simple")
  expect_true(!is.na(results))
  expect_length(results, 1)
  expect_gt(results, 0)
})
