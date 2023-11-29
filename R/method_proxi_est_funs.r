
#### Low-level estimation functions, written entirely by Andrew Ying
#### and taken (with permission) from his official github repository:
#### https://github.com/andrewyyp/Proximal_MSF
#### Only the formatting has been slightly changed by Robin Denz to fit
#### the rest of the package

censor_cumhaz_est <- function(censor_sorted_time, T0, Delta) {
  K_C <- length(censor_sorted_time)
  N <- length(T0)
  ## q value, N * 1
  IF <- matrix(0, nrow=K_C, ncol=N)
  for (j in seq_len(K_C)) {
    IF[j, ] <- (1 - Delta) * (T0 == censor_sorted_time[j]) /
      sum(T0 >= censor_sorted_time[j])
  }
  IF <- apply(IF, 2, cumsum)
  return(IF)
}

noncensor_cumhaz_compute <- function(sorted_time, censor_sorted_time, cum_haz) {
  noncensor_cumhaz <- c()
  K = length(sorted_time)
  K_C = length(censor_sorted_time)
  censor_sorted_time <- c(0, censor_sorted_time)
  cum_haz <- c(0, cum_haz)
  l = 1
  j = 1
  while(j <= K) {
    if(l + 1 <= K_C) {
      if(sorted_time[j] < censor_sorted_time[l + 1])
        noncensor_cumhaz <- c(noncensor_cumhaz, cum_haz[l])
      else {
        l = l + 1
        next
      }
    } else {
      noncensor_cumhaz <- c(noncensor_cumhaz, cum_haz[l])
    }
    j = j + 1
  }
  return(noncensor_cumhaz)
}

noncensor_cumhaz_IF_compute <- function(sorted_time, censor_sorted_time,
                                        cum_haz_IF) {
  noncensor_cumhaz_IF <- c()
  K <- length(sorted_time)
  K_C <- length(censor_sorted_time)
  censor_sorted_time <- c(0, censor_sorted_time)
  cum_haz_IF <- rbind(0, cum_haz_IF)
  l <- 1
  j <- 1
  while (j <= K) {
    if (l + 1 <= K_C) {
      if(sorted_time[j] < censor_sorted_time[l + 1]) {
        noncensor_cumhaz_IF <- rbind(noncensor_cumhaz_IF, cum_haz_IF[l, ])
      } else {
        l <- l + 1
        next
      }
    } else {
      noncensor_cumhaz_IF <- rbind(noncensor_cumhaz_IF, cum_haz_IF[l, ])
    }
    j <- j + 1
  }
  return(noncensor_cumhaz_IF)
}

hbridge <- function(b, h_W, h_Z, T0, sorted_time, noncensor_cumhaz) {
  K <- length(sorted_time)
  g <- 0

  for(j in seq_len(K)) {
    res <- cbind(h_Z, sorted_time[j]) / K
    integrand <- (T0 > sorted_time[j]) / exp(-noncensor_cumhaz[j]) -
      exp(-cbind(h_W, sorted_time[j]) %*% b * sorted_time[j])
    g <- g + c(integrand) * res
  }
  g <- colMeans(g)
  gmmf <- sum(g^2)
  return(gmmf)
}

hbridge_tmp <- function(b, h_W, h_Z, T0, sorted_time, noncensor_cumhaz) {
  K <- length(sorted_time)
  g <- 0

  for(j in seq_len(K)) {
    res <- cbind(h_Z, sorted_time[j]) / K
    integrand <- (T0 > sorted_time[j]) / exp(-noncensor_cumhaz[j]) -
      exp(-cbind(h_W, sorted_time[j]) %*% b * sorted_time[j])
    g <- g + c(integrand) * res
  }
  g <- colMeans(g)
  return(g)
}

h_IF <- function(b, h_W, h_Z, T0, sorted_time, noncensor_cumhaz,
                 noncensor_cumhaz_IF) {

  design_mat <- numDeriv::jacobian(
    func=function(...) {return(hbridge_tmp(...))},
    x=b, h_W=h_W, h_Z=h_Z, T0=T0,
    sorted_time=sorted_time,
    noncensor_cumhaz=noncensor_cumhaz
  )
  censor_mat <- h_censor(b, h_Z, T0, sorted_time, noncensor_cumhaz)
  K <- length(sorted_time)
  N <- nrow(h_W)
  eps <- 0

  for(j in seq_len(K)) {
    res <- cbind(h_Z, sorted_time[j]) / K
    integrand <- (T0 > sorted_time[j]) / exp(-noncensor_cumhaz[j]) -
      exp(-cbind(h_W, sorted_time[j]) %*% b * sorted_time[j])
    eps <- eps + c(integrand) * res
  }
  eps <- eps / N

  out <- -solve(design_mat) %*% (t(eps) + censor_mat %*% noncensor_cumhaz_IF)

  return(out)
}

h_censor <- function(b, h_Z, T0, sorted_time, noncensor_cumhaz) {
  K <- length(sorted_time)
  G <- matrix(0, nrow=length(b), ncol=K)
  for (j in seq_len(K)) {
    res <- cbind(h_Z, sorted_time[j]) / K
    integrand <- (T0 > sorted_time[j]) / exp(-noncensor_cumhaz[j])
    G[, j] <- colMeans(c(integrand) * res)
  }
  return(G)
}

MSE_func_q <- function(bridge_func, para, q_Y, q_W, q_Z){
  g0 <- bridge_func(para=para, q_Y=q_Y, q_W=q_W, q_Z=q_Z)
  g <- apply(g0, 2, mean)
  gmmf <- sum(g^2)
  return(gmmf)
}

qbridge <- function(para, q_Y, q_W, q_Z) {
  tlink <- 1 + exp(q_Z %*% para)
  g0 <- q_W
  g <- c(tlink) * g0 - q_Y
  return(g)
}

q_IF <- function(t, q_Y, q_W, q_Z) {
  out <- -solve(t(q_W) %*% (c(exp(q_Z %*% t)) * q_Z)) %*%
    t(c(1 + exp(q_Z %*% t)) * q_W - q_Y)
  return(out)
}

PIPW_surv <- function(sorted_time, time, noncensor_cumhaz, A, q_Z, t, a=1) {
  K <- length(sorted_time)
  N <- length(time)
  ## q value, N * 1
  q <- c(1 + exp(q_Z %*% t))
  IF <- matrix(0, nrow=K, ncol=N)
  for(j in seq_len(K)) {
    dN <- (A == a) * q * (time > sorted_time[j]) / exp(-noncensor_cumhaz[j])
    IF[j, ] <- dN / N
  }
  return(IF)
}

PIPW_q_dot <- function(sorted_time, time, noncensor_cumhaz, A, q_Z, t, a) {
  G <- numDeriv::jacobian(func=function(...) {return(rowSums(PIPW_surv(...)))},
                          x=t, sorted_time=sorted_time, time=time,
                          noncensor_cumhaz=noncensor_cumhaz,
                          A=A, q_Z=q_Z, a=a)
  return(G)
}

PIPW_censor_dot <- function(sorted_time, time, noncensor_cumhaz,
                            A, q_Z, t, a) {
  K <- length(sorted_time)
  N <- length(time)
  ## q value, N * 1
  q <- c(1 + exp(q_Z %*% t))
  G <- matrix(0, nrow=K, ncol=K)
  for(j in seq_len(K)) {
    dN <- (A == a) * q * (time > sorted_time[j]) / exp(-noncensor_cumhaz[j])
    G[j, j] <- mean(dN)
  }
  return(G)
}

PDR_surv <- function(b, t, time, sorted_time, noncensor_cumhaz,
                     h_W, q_Z, A, a=1) {
  K <- length(sorted_time)
  N <- nrow(h_W)
  IF <- matrix(0, nrow = K, ncol = N)
  q <- c(1 + exp(q_Z %*% t))

  for(j in seq_len(K)) {
    dh <- exp(-cbind(h_W, sorted_time[j]) %*% b * sorted_time[j])
    dN <- (A == a) * q * ((time > sorted_time[j]) /
                            exp(-noncensor_cumhaz[j]) - dh)
    IF[j, ] <- c(dN + dh)/N
  }
  return(IF)
}

PDR_h_dot <- function(b, t, time, sorted_time, noncensor_cumhaz,
                      h_W, q_Z, A, a) {

  G <- numDeriv::jacobian(func=function(...) {return(rowSums(PDR_surv(...)))},
                          x=b, sorted_time=sorted_time, time=time,
                          noncensor_cumhaz=noncensor_cumhaz,
                          A=A, h_W=h_W, t=t, q_Z=q_Z, a=a)
  return(G)
}

PDR_q_dot <- function(b, t, time, sorted_time, noncensor_cumhaz,
                      h_W, q_Z, A, a) {

  G <- numDeriv::jacobian(func=function(...) {return(rowSums(PDR_surv(...)))},
                          x=t, sorted_time=sorted_time, time=time,
                          noncensor_cumhaz=noncensor_cumhaz,
                          A=A, h_W=h_W, b=b, q_Z=q_Z, a=a)
  return(G)
}

PDR_censor_dot <- function(b, t, time, sorted_time, noncensor_cumhaz,
                           h_W, q_Z, A, a=1) {

  K <- length(sorted_time)
  N <- nrow(h_W)
  IF <- matrix(0, nrow=N, ncol=K)
  q <- c(1 + exp(q_Z %*% t))
  G <- matrix(0, nrow=K, ncol=K)
  for(j in seq_len(K)) {
    dh <- pmin(exp(-cbind(h_W, sorted_time[j]) %*% b * sorted_time[j]), 1)
    dN <- (A == a) * q * ((time > sorted_time[j]) /
                            exp(-noncensor_cumhaz[j]) - dh)
    G[j, j] <- mean(dN + dh)
  }
  return(G)
}

#### Additional high-level estimation functions written by Robin Denz
#### but also based on code from Andrew Ying below

## utility function to obtain the cumulative hazard from scratch
get_noncensor_cumhaz <- function(data, ev_time, event, sorted_time) {

  # (sorted) times at which individuals were censored
  censor_sorted_time <- sort(unique(data[, ev_time][data[, event]==0]))

  censor_cumhaz_sub <- censor_cumhaz_est(censor_sorted_time=censor_sorted_time,
                                         T0=data[, ev_time],
                                         Delta=data[, event])

  censor_cumhaz <- rowSums(censor_cumhaz_sub)
  censor_cumhaz_IF <- censor_cumhaz_sub - censor_cumhaz / nrow(data)

  noncensor_cumhaz <- noncensor_cumhaz_compute(
    sorted_time=sorted_time,
    censor_sorted_time=censor_sorted_time,
    cum_haz=censor_cumhaz
  )
  noncensor_cumhaz_IF <- noncensor_cumhaz_IF_compute(
    sorted_time=sorted_time,
    censor_sorted_time=censor_sorted_time,
    cum_haz_IF=censor_cumhaz_IF
  )

  # put important stuff together
  out <- list(noncensor_cumhaz=noncensor_cumhaz,
              noncensor_cumhaz_IF=noncensor_cumhaz_IF)

  return(out)
}

## utility function to fit q confounding bridge
get_q_confounding_bridge <- function(data, variable, adjust_vars,
                                     treatment_proxy, outcome_proxy,
                                     optim_method, optim_control) {

  # create design matrix to allow categorical adjust_vars
  all_vars <- c(variable, treatment_proxy, outcome_proxy, adjust_vars)
  x_dat <- as.data.frame(data[, all_vars])

  form <- paste0("~ ", paste0(all_vars, collapse=" + "))
  form <- remove_plus_at_end(form)

  mod_mat <- stats::model.matrix(stats::as.formula(form), data=x_dat)
  mod_mat <- as.matrix(mod_mat[,seq(2, ncol(mod_mat))])

  # number of parameters
  n_par <- ncol(mod_mat)

  # initial values for matrix generation / optimization
  inioptim_t <- ini_q <- rep(0, n_par)
  ini_q[2] <- 1

  # values for optimization
  q_Y <- matrix(rep(ini_q, each=nrow(data)), nrow=nrow(data), ncol=n_par)
  q_Z <- (-1)^(1 - data[, variable]) *
    cbind(1, mod_mat[, c(1, 2, 4:ncol(mod_mat))])
  q_W <- (-1)^(1 - data[, variable]) *
    cbind(1, mod_mat[, c(1, 3:ncol(mod_mat))])

  # call optim
  qlink <- stats::optim(par=inioptim_t, fn=MSE_func_q,
                        bridge_func=qbridge, q_Y=q_Y, q_W=q_W, q_Z=q_Z,
                        method=optim_method, hessian=FALSE,
                        control=optim_control)

  # get influence function
  t_IF <- q_IF(t=qlink$par, q_Y=q_Y, q_W=q_W, q_Z=q_Z)

  # put together
  out <- list(q_Z=q_Z,
              qlink=qlink,
              t=qlink$par,
              t_IF=t_IF)

  return(out)
}

## utility function to obtain the h confounding bridge
get_h_confounding_bridge <- function(data, variable, adjust_vars,
                                     treatment_proxy, outcome_proxy,
                                     noncens, sorted_time, ev_time,
                                     optim_method, optim_control) {

  # create design matrix to allow categorical adjust_vars
  all_vars <- c(variable, treatment_proxy, outcome_proxy, adjust_vars)
  x_dat <- as.data.frame(data[, all_vars])

  form <- paste0("~ ", paste0(all_vars, collapse=" + "))
  form <- remove_plus_at_end(form)

  mod_mat <- stats::model.matrix(stats::as.formula(form), data=x_dat)
  mod_mat <- as.matrix(mod_mat[,seq(2, ncol(mod_mat))])

  # initialize starting values
  n_par <- ncol(mod_mat) + 1
  inioptim_b <- rep(0.1, n_par)

  h_W_1 <- as.matrix(cbind(1, 1, mod_mat[, 3:ncol(mod_mat)]))
  h_W_0 <- as.matrix(cbind(1, 0, mod_mat[, 3:ncol(mod_mat)]))

  h_W <- as.matrix(cbind(1, mod_mat[, c(1, 3:ncol(mod_mat))]))
  h_Z <- as.matrix(cbind(1, mod_mat[, c(1, 2, 4:ncol(mod_mat))]))

  # optimize them
  hlink <- stats::optim(par=inioptim_b, fn=hbridge, h_W=h_W, h_Z=h_Z,
                        T0=data[, ev_time], sorted_time=sorted_time,
                        noncensor_cumhaz=noncens$noncensor_cumhaz,
                        method=optim_method, hessian=FALSE,
                        control=optim_control)

  b <- hlink$par

  # extract the IF of q function
  b_IF <- h_IF(b=b, h_W=h_W, h_Z=h_Z, T0=data[, ev_time],
               sorted_time=sorted_time,
               noncensor_cumhaz=noncens$noncensor_cumhaz,
               noncensor_cumhaz_IF=noncens$noncensor_cumhaz_IF)

  # put together
  out <- list(h_Z=h_Z,
              h_W=h_W,
              h_W_1=h_W_1,
              h_W_0=h_W_0,
              hlink=hlink,
              b=hlink$par,
              b_IF=b_IF)

  return(out)
}
