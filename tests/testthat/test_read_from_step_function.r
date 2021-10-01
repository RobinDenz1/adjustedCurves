
test_dat <- data.frame(time=c(0, 0.1, 0.2, 0.3, 0.6, 0.63),
                       surv=c(1, 0.99, 0.8, 0.76, 0.5, 0.3),
                       cif=1-c(1, 0.99, 0.8, 0.76, 0.5, 0.3))

test_that("surv, t = 0, start", {
  expect_equal(adjustedCurves:::read_from_step_function(0, step_data=test_dat,
                                                        est="surv"),
               1)
})

test_that("cif, t = 0, start", {
  expect_equal(adjustedCurves:::read_from_step_function(0, step_data=test_dat,
                                                        est="cif"),
               0)
})

test_that("surv, t = 0.2, on jump", {
  expect_equal(adjustedCurves:::read_from_step_function(0.2, step_data=test_dat,
                                                        est="surv"),
               0.8)
})

test_that("cif, t = 0.2, on jump", {
  expect_equal(adjustedCurves:::read_from_step_function(0.2, step_data=test_dat,
                                                        est="cif"),
               0.2)
})

test_that("surv, t = 0.4, in slope = 0", {
  expect_equal(adjustedCurves:::read_from_step_function(0.4, step_data=test_dat,
                                                        est="surv"),
               0.76)
})

test_that("cif, t = 0.4, in slope = 0", {
  expect_equal(adjustedCurves:::read_from_step_function(0.4, step_data=test_dat,
                                                        est="cif"),
               0.24)
})


test_that("surv, t = 0.63, end", {
  expect_equal(adjustedCurves:::read_from_step_function(0.63, step_data=test_dat,
                                                        est="surv"),
               0.3)
})

test_that("cif, t = 0.63, end", {
  expect_equal(adjustedCurves:::read_from_step_function(0.63, step_data=test_dat,
                                                        est="cif"),
               0.7)
})

test_that("surv, t = 0.64, no extrapolation", {
  expect_true(is.na(adjustedCurves:::read_from_step_function(0.64,
                                                             step_data=test_dat,
                                                             est="surv")))
})

test_that("cif, t = 0.64, no extrapolation", {
  expect_true(is.na(adjustedCurves:::read_from_step_function(0.64,
                                                             step_data=test_dat,
                                                             est="cif")))
})
