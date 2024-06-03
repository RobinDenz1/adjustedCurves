
para_set <- list(mu_X=1.1,
                 sigma_X=0.75,
                 mu_U=1.1,
                 sigma_U=0.75,
                 alpha_A=c(0.3, 0.4, -0.6),
                 mu_Z=c(-0.2, -0.3, 0.65),
                 sigma_Z=0.5,
                 mu_W=c(-0.6, 0.4, 0.65),
                 sigma_W=0.5,
                 mu_T0=c(0.1, 0.6, 0.25, 0.5),
                 mu_C=0.2,
                 admin_C=2
)

data_gen <- function(N, para_set, a = NULL) {
  # generate X, U
  X <- para_set$mu_X + rnorm(N, 0, para_set$sigma_X)
  U <- para_set$mu_U + rnorm(N, 0, para_set$sigma_U)
  X <- pmax(X, 0)
  U <- pmax(U, 0)

  if (is.null(a)) {
    # generate A
    prop_score_0 <- 1/(1 + exp(-cbind(1, X, U) %*% para_set$alpha_A))
    A <- rbinom(N, 1, prop_score_0)
  } else {
    A <- rep(a, N)
  }
  # generate Z
  Z <- cbind(1, X, U) %*% para_set$mu_Z + rnorm(N, 0, para_set$sigma_Z)

  # generate W
  W <- cbind(1, X, U) %*% para_set$mu_W + rnorm(N, 0, para_set$sigma_W)

  # generate Y
  T0 <- rexp(N, rate = cbind(1, A, X, U) %*% para_set$mu_T0)

  C <- rexp(N, rate = para_set$mu_C)
  C <- pmin(C, para_set$admin_C)
  if (is.null(a)) {
    df <- data.frame(X, U, A, Z, W, T0 = pmin(T0, C), Delta = (T0 <= C))
  } else {
    df <- data.frame(X, U, A, Z, W, T0 = T0, Delta = rep(1, N))
  }
  return(df)
}

set.seed(4356)
n_sim <- 100
dat <- data_gen(n_sim, para_set)
dat$A <- factor(dat$A)
dat$X_2 <- rnorm(n_sim)
dat$X_3 <- factor(sample.int(3, size=n_sim, replace=TRUE))

test_that("default arguments", {
  out <- adjustedsurv(data=dat,
                      variable="A",
                      ev_time="T0",
                      event="Delta",
                      times=NULL,
                      adjust_vars="X",
                      treatment_proxy="Z",
                      outcome_proxy="W",
                      conf_int=FALSE,
                      method="prox_aiptw")
  expect_s3_class(out, "adjustedsurv")
  expect_true(is.numeric(out$adj$surv))
})

test_that("using conf_int", {
  out <- adjustedsurv(data=dat,
                      variable="A",
                      ev_time="T0",
                      event="Delta",
                      times=NULL,
                      adjust_vars="X",
                      treatment_proxy="Z",
                      outcome_proxy="W",
                      conf_int=TRUE,
                      method="prox_aiptw")
  expect_s3_class(out, "adjustedsurv")
  expect_true(is.numeric(out$adj$surv))
  expect_true(is.numeric(out$adj$ci_lower))
})

test_that("using bootstrapping", {
  out <- adjustedsurv(data=dat,
                      variable="A",
                      ev_time="T0",
                      event="Delta",
                      times=NULL,
                      adjust_vars="X",
                      treatment_proxy="Z",
                      outcome_proxy="W",
                      conf_int=FALSE,
                      method="prox_aiptw",
                      bootstrap=TRUE,
                      n_boot=2)
  expect_s3_class(out, "adjustedsurv")
  expect_true(is.numeric(out$adj$surv))
  expect_true(is.numeric(out$boot_adj$ci_lower))
})

test_that("with times", {
  out <- adjustedsurv(data=dat,
                      variable="A",
                      ev_time="T0",
                      event="Delta",
                      times=c(0.5, 1),
                      adjust_vars="X",
                      treatment_proxy="Z",
                      outcome_proxy="W",
                      conf_int=FALSE,
                      method="prox_aiptw")
  expect_s3_class(out, "adjustedsurv")
  expect_true(is.numeric(out$adj$surv))
  expect_true(length(out$adj$time)==4)
})

test_that("changing some custom arguments", {
  out <- adjustedsurv(data=dat,
                      variable="A",
                      ev_time="T0",
                      event="Delta",
                      times=NULL,
                      adjust_vars="X",
                      treatment_proxy="Z",
                      outcome_proxy="W",
                      conf_int=FALSE,
                      method="prox_aiptw",
                      optim_method_q="Nelder-Mead",
                      optim_control_q=list(maxit=40),
                      return_fit=FALSE)
  expect_s3_class(out, "adjustedsurv")
  expect_true(is.numeric(out$adj$surv))
})

test_that("with multiple numeric adjust_vars", {
  out <- adjustedsurv(data=dat,
                      variable="A",
                      ev_time="T0",
                      event="Delta",
                      adjust_vars=c("X", "X_2"),
                      treatment_proxy="Z",
                      outcome_proxy="W",
                      conf_int=FALSE,
                      method="prox_aiptw")
  expect_s3_class(out, "adjustedsurv")
  expect_true(is.numeric(out$adj$surv))
})

test_that("with multiple partially factor adjust_vars", {
  out <- adjustedsurv(data=dat,
                      variable="A",
                      ev_time="T0",
                      event="Delta",
                      adjust_vars=c("X", "X_2", "X_3"),
                      treatment_proxy="Z",
                      outcome_proxy="W",
                      conf_int=FALSE,
                      method="prox_aiptw")
  expect_s3_class(out, "adjustedsurv")
  expect_true(is.numeric(out$adj$surv))
})
