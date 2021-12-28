
weights <- c(0, 1, 0.3, 0.4, 19)

trimmed <- trim_weights(weights, trim=18)

test_that("trimming", {
  expect_lt(max(trimmed), 19)
})
