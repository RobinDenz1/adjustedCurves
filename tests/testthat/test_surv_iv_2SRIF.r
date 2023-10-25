
set.seed(42)

## This example has been (mostly) taken from the online supplement of the
## original article by Martinez-Camblor et al. (2021)

# generate some data
n <- 1000
t <- seq(0, 10, 0.01)
bu <- log(2)
hr <- 2
v <- 2
a <- 1

U <- stats::rnorm(n)
Z <- stats::rnorm(n)
Z_2 <- stats::rnorm(n)
Z_3 <- factor(sample(c("A", "B", "C"), size=n, replace=TRUE))
W <- stats::rnorm(n)
e <- stats::rnorm(n)

X0 <- (U + Z + a*W + (v - a^2)^0.5*e >= 0)
L0 <- Z + bu*U
L1 <- log(hr) + Z + bu*U
T <- stats::rexp(n, 0.005)

T0 <- T/exp(L0)
T1 <- T/exp(L1)

censor <- stats::rexp(n, 0.05)
time1 <- pmin(ifelse(X0==1,T1,T0), censor)
status1 <- 1-(censor==time1)
time <- pmin(time1, 10)
status <- ifelse(time1 > 10, 0, status1)

dt <- as.data.frame(cbind(time, status, X0, Z, Z_2, Z_3, W))
dt$X0 <- factor(dt$X0)
dt$Z_3 <- factor(dt$Z_3)

test_that("general case", {
  adjsurv <- adjustedsurv(data=dt,
                          variable="X0",
                          ev_time="time",
                          event="status",
                          method="iv_2SRIF",
                          adjust_vars="Z",
                          instrument="W")
  expect_s3_class(adjsurv, "adjustedsurv")
  expect_true(is.numeric(adjsurv$adjsurv$surv))
})

test_that("with bootstrapping", {
  adjsurv <- adjustedsurv(data=dt,
                          variable="X0",
                          ev_time="time",
                          event="status",
                          method="iv_2SRIF",
                          adjust_vars="Z",
                          instrument="W",
                          bootstrap=TRUE,
                          n_boot=2)
  expect_s3_class(adjsurv, "adjustedsurv")
  expect_true(is.numeric(adjsurv$adjsurv$surv))
})

test_that("multiple adjust_vars", {
  adjsurv <- adjustedsurv(data=dt,
                          variable="X0",
                          ev_time="time",
                          event="status",
                          method="iv_2SRIF",
                          adjust_vars=c("Z", "Z_2"),
                          instrument="W")
  expect_s3_class(adjsurv, "adjustedsurv")
  expect_true(is.numeric(adjsurv$adjsurv$surv))
})

test_that("multiple adjust_vars with a factor variable", {
  adjsurv <- adjustedsurv(data=dt,
                          variable="X0",
                          ev_time="time",
                          event="status",
                          method="iv_2SRIF",
                          adjust_vars=c("Z", "Z_2", "Z_3"),
                          instrument="W")
  expect_s3_class(adjsurv, "adjustedsurv")
  expect_true(is.numeric(adjsurv$adjsurv$surv))
})

test_that("no adjust_vars", {
  adjsurv <- adjustedsurv(data=dt,
                          variable="X0",
                          ev_time="time",
                          event="status",
                          method="iv_2SRIF",
                          adjust_vars=NULL,
                          instrument="W")
  expect_s3_class(adjsurv, "adjustedsurv")
  expect_true(is.numeric(adjsurv$adjsurv$surv))
})
