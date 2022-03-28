
test_dat <- data.frame(time=c(0, 0.1, 0.2, 0.3, 0.6,
                              0, 0.15, 0.25, 0.35, 0.7),
                       surv=c(1, 0.9, 0.8, 0.75, 0.5,
                              1, 0.8, 0.7, 0.65, 0.3),
                       group=c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1))
test_dat$group <- factor(test_dat$group, levels=c("1", "0"))

results <- data.frame(time=c(0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.6, 0.7),
                      surv=c(0, 0.1, -0.1, 0, -0.1, -0.05, -0.1, 0.15, NA))


test_that("example test, same length", {
  out <- difference_function(test_dat, times=sort(unique(test_dat$time)))
  expect_equal(out, results)
})

test_dat <- data.frame(time=c(0, 0.1, 0.2, 0.3, 0.6,
                              0, 0.15, 0.25, 0.35, 0.7, 0.71),
                       surv=c(1, 0.9, 0.8, 0.75, 0.5,
                              1, 0.8, 0.7, 0.65, 0.3, 0.29),
                       group=c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1))
test_dat$group <- factor(test_dat$group, levels=c("1", "0"))

results <- data.frame(time=c(0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.6, 0.7,
                             0.71),
                      surv=c(0, 0.1, -0.1, 0, -0.1, -0.05, -0.1, 0.15, NA, NA))

test_that("example test, different length", {
  out <- difference_function(test_dat, times=sort(unique(test_dat$time)))
  expect_equal(out, results)
})
