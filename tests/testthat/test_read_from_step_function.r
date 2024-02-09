
test_dat <- data.frame(time=c(0, 0.1, 0.2, 0.3, 0.6, 0.63),
                       surv=c(1, 0.99, 0.8, 0.76, 0.5, 0.3),
                       cif=1-c(1, 0.99, 0.8, 0.76, 0.5, 0.3))

test_that("surv, t = 0, start", {
  out <- read_from_step_function(0, data=test_dat, est="surv")
  expect_equal(out, 1)
})

test_that("cif, t = 0, start", {
  out <- read_from_step_function(0, data=test_dat, est="cif")
  expect_equal(out, 0)
})

test_that("surv, t = 0.2, on jump", {
  out <- read_from_step_function(0.2, data=test_dat, est="surv")
  expect_equal(out, 0.8)
})

test_that("cif, t = 0.2, on jump", {
  out <- read_from_step_function(0.2, data=test_dat, est="cif")
  expect_equal(out, 0.2)
})

test_that("surv, t = 0.4, in slope = 0", {
  out <- read_from_step_function(0.4, data=test_dat, est="surv")
  expect_equal(out, 0.76)
})

test_that("cif, t = 0.4, in slope = 0", {
  out <- read_from_step_function(0.4, data=test_dat, est="cif")
  expect_equal(out, 0.24)
})

test_that("surv, t = 0.63, end", {
  out <- read_from_step_function(0.63, data=test_dat, est="surv")
  expect_equal(out, 0.3)
})

test_that("cif, t = 0.63, end", {
  out <- read_from_step_function(0.63, data=test_dat, est="cif")
  expect_equal(out, 0.7)
})

test_that("surv, t = 0.64, no extrapolation", {
  out <- read_from_step_function(0.64, data=test_dat, est="surv")
  expect_true(is.na(out))
})

test_that("cif, t = 0.64, no extrapolation", {
  out <- read_from_step_function(0.64, data=test_dat, est="cif")
  expect_true(is.na(out))
})

test_dat <- data.frame(time=c(0, 0.1, 0.2, 0.3, 0.6, 0.63),
                       surv=c(1, 0.99, 0.8, 0.76, 0.5, 0),
                       cif=1-c(1, 0.99, 0.8, 0.76, 0.5, 0))

test_that("surv, t = 0.64, valid extrapolation", {
  out <- read_from_step_function(1.64, data=test_dat, est="surv")
  expect_equal(out, 0)
})

test_that("cif, t = 0.64, valid extrapolation", {
  out <- read_from_step_function(1.64, data=test_dat, est="cif")
  expect_equal(out, 1)
})
