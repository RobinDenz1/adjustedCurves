
### Step Function

test_dat_surv <- data.frame(time=c(0, 0.1, 0.2, 0.3, 0.6),
                            surv=c(1, 0.9, 0.8, 0.75, 0.5))
test_dat_cif <- data.frame(time=c(0, 0.1, 0.2, 0.3, 0.6),
                           cif=1-c(1, 0.9, 0.8, 0.75, 0.5))

test_that("surv, full", {
  out <- exact_integral(test_dat_surv,
                        from=0,
                        to=0.6,
                        est="surv",
                        interpolation="steps")
  expect_equal(out, 0.495)
})

test_that("cif, full", {
  out <- exact_integral(test_dat_cif,
                        from=0,
                        to=0.6,
                        est="cif",
                        interpolation="steps")
  expect_equal(out, 0.105)
})

test_that("surv, to = 0.2", {
  out <- exact_integral(test_dat_surv,
                        from=0,
                        to=0.2,
                        est="surv",
                        interpolation="steps")
  expect_equal(out, 0.19)
})

test_that("surv, to = c(0.2, 0.4)", {
  out <- exact_integral(test_dat_surv,
                        from=0,
                        to=c(0.2, 0.4),
                        est="surv",
                        interpolation="steps")
  expect_equal(out, c(0.19, 0.345))
})

test_that("cif, to = 0.2", {
  out <- exact_integral(test_dat_cif,
                        from=0,
                        to=0.2,
                        est="cif",
                        interpolation="steps")
  expect_equal(out, 0.01)
})

test_that("cif, to = c(0.2, 0.4)", {
  out <- exact_integral(test_dat_cif,
                        from=0,
                        to=c(0.2, 0.4),
                        est="cif",
                        interpolation="steps")
  expect_equal(out, c(0.01, 0.055))
})

test_that("surv, from = 0.2", {
  out <- exact_integral(test_dat_surv,
                        from=0.2,
                        to=0.6,
                        est="surv",
                        interpolation="steps")
  expect_equal(out, 0.305)
})

test_that("cif, from = 0.2", {
  out <- exact_integral(test_dat_cif,
                        from=0.2,
                        to=0.6,
                        est="cif",
                        interpolation="steps")
  expect_equal(out, 0.095)
})

test_that("surv, from = 0.2, to = 0.4", {
  out <- exact_integral(test_dat_surv,
                        from=0.2,
                        to=0.4,
                        est="surv",
                        interpolation="steps")
  expect_equal(out, 0.155)
})

test_that("surv, from = 0.2, to = c(0.4, 0.6)", {
  out <- exact_integral(test_dat_surv,
                        from=0.2,
                        to=c(0.4, 0.6),
                        est="surv",
                        interpolation="steps")
  expect_equal(out, c(0.155, 0.305))
})

test_that("cif, from = 0.2, to = 0.4", {
  out <- exact_integral(test_dat_cif,
                        from=0.2,
                        to=0.4,
                        est="cif",
                        interpolation="steps")
  expect_equal(out, 0.045)
})

test_that("cif, from = 0.2, to = c(0.4, 0.6)", {
  out <- exact_integral(test_dat_cif,
                        from=0.2,
                        to=c(0.4, 0.6),
                        est="cif",
                        interpolation="steps")
  expect_equal(out, c(0.045, 0.095))
})
