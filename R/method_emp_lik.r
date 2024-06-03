# Copyright (C) 2021  Robin Denz
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

## All functions in this file were written by Fangfang Bai and Xiaofei Wang
## and were taken directly from the official github repository
## See: https://github.com/kimihua1995
## The only exception is the "surv_emp_lik" function.
## Since that package is published under a General Public License, this is
## perfectly legal. Permission from the authors was also granted.

# CHANGES BY ROBIN DENZ (21.05.2021):
#   - minor changes to the code layout
#   - changed 1:n in for loop to seq_len(n)
jacob.fun <- function(treat, x, psix, a, tau, d) {

  n <- length(treat)
  n1 <- sum(treat)
  n2 <- n - n1
  psi <- as.matrix(psix - t(array(rep(a, n), dim=c(d, n))))
  midd1 <- (n1 + psi %*% tau)^2
  midd11 <- (n1 + psi %*% tau)
  midd2 <- (n2 - psi %*% tau)^2
  midd22 <- (n2 - psi %*% tau)
  E <- diag(rep(1, d))
  a1 <- a2 <- a3 <- a4 <-0
  for (i in seq_len(n)) {
    di <- treat[i]
    psii <- as.matrix(psix[i,] - a)
    midd1i <- midd1[i]
    midd12i <- midd11[i]
    midd2i <- midd2[i]
    midd23i <- midd22[i]
    a1 <- a1 + (-E * midd12i + psii %*% t(tau)) * di / midd1i
    a2 <- a2 + (-psii %*% t(psii)) * di / midd1i
    a3 <- a3 + (-E * midd23i - psii %*% t(tau)) * (1 - di) / midd2i
    a4 <- a4 + (psii %*% t(psii)) * (1 - di) / midd2i
  }
  a1 <- a1
  a2 <- a2
  a3 <- a3
  a4 <- a4
  Jacob <- rbind(cbind(a1, a2), cbind(a3, a4))

}

## This is the function for estimating a and tau in EL likelihood function
# CHANGES BY ROBIN DENZ (21.05.2021):
#   - minor changes to the code layout
est.a.tau <- function(treat, x, psix, a, tau, d) {

  n <- length(treat)
  n1 <- sum(treat)
  n2 <- n - n1
  psi <- psix - t(array(rep(a, n), dim=c(d, n)))
  rep_treat <- array(rep(treat, d), dim=c(n, d))
  f1 <- apply(rep_treat * psi / matrix(rep((n1 + psi %*% tau), d), ncol=d),
              2, sum)
  f2 <- apply((1 - rep_treat) * psi / matrix(rep((n2 - psi %*% tau), d),
                                             ncol=d),
              2, sum)
  f <- matrix(c(f1, f2), ncol=1)

}

## This is the estimating function to obtain the estimator
## of a and tau by netwon raphson method

# CHANGES BY ROBIN DENZ (21.05.2021):
#   - added newton_tol argument, making it possible for the user
#     to set a value for the newton-raphson tolerance
#   - added max_iter argument allowing the user to change the
#     maximum number of iterations allowed
#   - changed 1:max_iter in for loop to seq_len(max_iter)
#   - minor changes to the code layout
# CHANGES BY ROBIN DENZ (27.12.2021):
#   - added warning if the algorithm did not converge
#   - made sure this function still returns something in that case
estamtor <- function(treat, x, psix, a, tau, d, max_iter, newton_tol) {

  n <- length(treat)
  n1 <- sum(treat)
  n2 <- n-n1
  para <- rbind(a, tau)
  a0 <- as.matrix(para[1:d])
  tau0 <- as.matrix(para[-(1:d)])
  fvalue <- est.a.tau(treat, x, psix, a0, tau0, d)

  converge <- FALSE
  if(max(abs(fvalue)) < newton_tol) {

    converge <- TRUE
    return(list(a=as.matrix(para[1:d]),
                tau=as.matrix(para[-(1:d)]),
                converge=converge))
  }

  dfvalue <- jacob.fun(treat, x, psix, a0, tau0, d)

  for (i in seq_len(max_iter)) {

    para1 <- para - MASS::ginv(dfvalue) %*% fvalue

    error <- max(abs(para1 - para))
    if (error < newton_tol | max(abs(fvalue)) < newton_tol) {
      converge <- TRUE
      return(list(a=as.matrix(para1[1:d]),
                  tau=as.matrix(para1[-(1:d)]),
                  converge=converge))
    }

    para <- para1
    a0 <- as.matrix(para[1:d])
    tau0 <- as.matrix(para[-(1:d)])
    fvalue <- est.a.tau(treat, x, psix, a0, tau0, d)
    dfvalue <- jacob.fun(treat, x, psix, a0, tau0, d)
  }

  warning("Algorithm did not converge.",
          " Increasing the max_iter value or decreasing the newton_tol",
          " value might help.",
          call.=FALSE)

  return(list(a=as.matrix(para1[1:d]),
              tau=as.matrix(para1[-(1:d)]),
              converge=converge))
}

# CHANGES ROBIN DENZ (21.05.2021):
#   - added max_iter, newton_tol, see above
#   - minor changes to the code layout
estimator.pi <- function(y, delta, treat, x, psix, max_iter, newton_tol) {

  dd <- ncol(psix)
  tau0 <- array(rep(0, dd), dim=c(dd, 1))
  a0 <- array(rep(0, dd), dim=c(dd, 1))

  n <- length(y)
  n1 <- sum(treat)
  n2 <- n - n1
  # the estimator for a and tau
  estnusi <- estamtor(treat, x, psix, a0, tau0, dd, max_iter, newton_tol)
  pi <- rep(1, n)

  if (estnusi$converge==TRUE && length(estnusi$converge)!=0) {
    a <- estnusi$a
    tau <- estnusi$tau
    psi <- psix - t(array(rep(a, n), dim=c(dd, n)))
    pi1 <- 1 / (n1 + psi %*% tau)
    pi2 <- 1 / (n2 - psi %*% tau)
    pi[treat==1] <- pi1[treat==1]
    pi[treat==0] <- pi2[treat==0]
  }
  return(pi)
}

## main function
# CHANGES BY ROBIN DENZ (21.05.2021):
#   - removed get_sd argument and functionality because it is just
#     bootstrapping, which is already implemented in "adjustedsurv"
#   - added "gtol" argument, making it possible to choose a value
#     for the numeric tolerance
#   - added the line "pi[pi<0] <- gtol"
#   - minor changes to the code layout
#   - added max_iter, newton_tol, see above
#   - moved "require(MASS)" to load_needed_packages function, used
#     requireNamespace("MASS") instead
el.est <- function(y, delta, treat, x, psix_moment=c("first", "second"),
                   treat.select, t, standardize=FALSE, gtol=0.00001,
                   max_iter=100, newton_tol=1.0e-06){

  if (standardize) {
    x <- scale(x)
  }

  if (psix_moment=="first") {
    psix <- x
  } else if (psix_moment=="second") {
    p <- ncol(x)
    combb <- utils::combn(p, 2)
    res <- apply(x, 1, function(x) {apply(combb, 2, function(y) prod(x[y]))})
    psix <- cbind(x^2, t(res))
  }

  pi <- estimator.pi(y, delta, treat, x, psix, max_iter, newton_tol)
  # New line added to keep method from producing nonsense estimates
  # in some situations
  pi[pi<0] <- gtol

  y <- y[treat==treat.select]
  delta <- delta[treat==treat.select]
  pi <- pi[treat==treat.select]

  n <- length(y)
  m <- length(t)
  z0 <- y[delta==1]
  pi.1 <- pi[delta==1]

  N <- length(z0)
  if (N > 0) {
    II0 <- (matrix(rep(z0, m), ncol=m) <= matrix(rep(t, N), ncol=m,
                                                 byrow=TRUE)) *
      matrix(rep(pi.1, m), ncol=m)
    I0 <- as.numeric(matrix(rep(y, N), ncol=N) >=
                     matrix(rep(z0, n), ncol=N, byrow=TRUE)) *
      matrix(rep(pi, N), ncol=N)
    ssum <- matrix(rep(apply(I0, 2,sum), m), ncol=m)
    temp <- 1 - II0 / (ssum + gtol)
  } else {
    temp <- matrix(1, n, m)
  }

  suvdf <- apply(temp, 2, prod)

  return(suvdf)
}

## Using Empirical Likelihood Estimation
#' @export
surv_emp_lik <- function(data, variable, ev_time, event, conf_int=FALSE,
                         conf_level=0.95, times, treatment_vars, moment="first",
                         standardize=FALSE, gtol=0.00001,
                         max_iter=100, newton_tol=1.0e-06) {

  # turn it into numeric
  levs <- levels(data[, variable])
  data[, variable] <- ifelse(data[, variable]==levs[1], 0, 1)

  # create design matrix for function call
  x_dat <- as.data.frame(data[,treatment_vars])
  if (length(treatment_vars==1)) {
    colnames(x_dat) <- treatment_vars
  }
  form <- stats::as.formula(paste0("~ ",
                                   paste0(treatment_vars, collapse=" + ")))
  mod_mat <- stats::model.matrix(form, data=x_dat)
  mod_mat <- as.matrix(mod_mat[,seq(2, ncol(mod_mat))])

  # call function twice
  el_0 <- el.est(y=data[, ev_time],
                 delta=data[, event],
                 treat=data[, variable],
                 x=mod_mat,
                 treat.select=0,
                 t=times,
                 psix_moment=moment,
                 standardize=standardize,
                 gtol=gtol,
                 max_iter=max_iter,
                 newton_tol=newton_tol)

  el_1 <- el.est(y=data[, ev_time],
                 delta=data[, event],
                 treat=data[, variable],
                 x=mod_mat,
                 treat.select=1,
                 t=times,
                 psix_moment=moment,
                 standardize=standardize,
                 gtol=gtol,
                 max_iter=max_iter,
                 newton_tol=newton_tol)

  plotdata <- data.frame(time=c(times, times),
                         surv=c(el_0, el_1),
                         group=c(rep(levs[1], length(el_0)),
                                 rep(levs[2], length(el_0))))
  plotdata$group <- factor(plotdata$group, levels=levs)

  output <- list(plotdata=plotdata)
  class(output) <- "adjustedsurv.method"

  return(output)
}
