
test_that("all equal p-values should output that value", {
  expect_equal(pool_p_values(c(0, 0, 0, 0, 0)), 0)
})

test_that("all equal p-values should be invariant to tol", {
  expect_equal(pool_p_values(c(0, 0, 0, 0, 0), tol=0.1), 0)
})

test_that("no error when using different ones", {
  expect_error(pool_p_values(c(0.1, 0.14, 0.18, 0.07)), NA)
})

# very small, obviously not very valid simulation
set.seed(42)
pooled <- numeric(1000)
for (i in seq_len(1000)) {
  pooled[i] <- pool_p_values(p_values=stats::runif(n=10, min=0.4, max=0.6))
}

# results however should be very close to 0.5, never exceed 0 or 1
test_that("p-values should be <= 1", {
  expect_true(all(pooled <= 1))
})

test_that("p-values should be >= 0", {
  expect_true(all(pooled >= 0))
})

test_that("p-values should be close to 0.5 here", {
  expect_true(mean(pooled) > 0.49 & mean(pooled) < 0.51)
})
