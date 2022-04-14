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

## takes a value x at which to read from the step function
## and step function data from which to read it
read_from_step_function <- function(x, data, est="surv", time="time") {

  # keep only data with non-missing est
  data <- data[which(!is.na(data[, est])), ]

  # no extrapolation
  if (x > max(data[, time])) {
    return(NA)
  }

  # otherwise get value
  check <- data[which(data[, time] <= x), ]
  if (nrow(check)==0) {
    if (est=="surv") {
      val <- 1
    } else if (est=="cif") {
      val <- 0
    } else {
      val <- NA
    }
  } else {
    val <- check[, est][which(check[, time]==max(check[, time]))][1]
  }
  return(val)
}

## takes a value x at which to read from the linear function
## and step function data from which to read it
read_from_linear_function <- function(x, data, est="surv", time="time") {

  # keep only data with non-missing est
  data <- data[which(!is.na(data[, est])), ]

  time_vec <- data[, time]
  est_vec <- data[, est]

  # add zero to beginning if necessary
  if (est=="surv") {
    val <- suppressWarnings(stats::approx(x=time_vec, y=est_vec, xout=x,
                                          yleft=1)$y)
  } else if (est=="cif") {
    val <- suppressWarnings(stats::approx(x=time_vec, y=est_vec, xout=x,
                                          yleft=0)$y)
  } else {
    val <- suppressWarnings(stats::approx(x=time_vec, y=est_vec, xout=x)$y)
  }
  return(val)
}

## shortcut to one of the two functions above
read_from_fun <- function(x, data, interpolation, est="surv", time="time") {

  if (interpolation=="steps") {
    val <- read_from_step_function(x=x, data=data, est=est, time=time)
  } else if (interpolation=="linear") {
    val <- read_from_linear_function(x=x, data=data, est=est, time=time)
  }
  return(val)
}

## calculate integral of a linear function using the trapezoid rule
trapezoid_integral <- function(x, y) {
  area <- 0
  for (i in seq_len((length(x)-1))) {
    a <- x[i]
    b <- x[i+1]
    f_a <- y[i]
    f_b <- y[i+1]

    h <- b - a

    area <- area + ((h / 2) * (f_a + f_b))
  }
  return(area)
}

## calculate exact integral of a step function
stepfun_integral <- function(x, y) {
  area <- 0
  for (i in seq_len((length(x)-1))) {
    x1 <- x[i]
    x2 <- x[i+1]
    rect_area <- (x2 - x1) * y[i]
    area <- area + rect_area
  }
  return(area)
}

## calculate exact integral of either step or linear functions
## in a given interval
# NOTE: 'data' needs to be a data.frame (!) with columns 'time' and 'est',
# sorted by time with no duplicates in time
exact_integral <- function(data, from, to, est, interpolation) {

  # constrain function end
  latest <- read_from_fun(to, data=data, est=est, interpolation=interpolation)
  data <- data[data$time <= to, ]

  if (!to %in% data$time) {
    temp <- data.frame(time=to)
    temp[, est] <- latest
    data <- rbind(data, temp)
  }

  # constrain function beginning
  if (from != 0) {
    earliest <- read_from_fun(from, data=data, est=est,
                              interpolation=interpolation)
    data <- data[data$time >= from, ]

    if (!from %in% data$time) {
      temp <- data.frame(time=from)
      temp[, est] <- earliest
      data <- rbind(temp, data)
    }
  }

  if (anyNA(data)) {
    return(NA)
  }

  if (interpolation=="steps") {
    area <- stepfun_integral(x=data$time, y=data[, est])
  } else {
    area <- trapezoid_integral(x=data$time, y=data[, est])
  }
  return(area)
}

## calculate difference between two functions at arbitrary time values
## using either linear or step-function interpolation
difference_function <- function(adj, times, est="surv", interpolation="steps",
                                conf_int=FALSE, conf_level=0.95) {

  levs <- levels(adj$group)
  adjsurv_0 <- adj[which(adj$group==levs[1]), ]
  adjsurv_1 <- adj[which(adj$group==levs[2]), ]

  if (nrow(adjsurv_0)==nrow(adjsurv_1) && all(adjsurv_0$time==adjsurv_1$time) &&
      interpolation=="steps") {
    surv_0 <- adjsurv_0[, est]
    surv_1 <- adjsurv_1[, est]

    if (conf_int) {
      se_0 <- adjsurv_0$se
      se_1 <- adjsurv_1$se
    }
  } else {
    surv_0 <- vapply(times, read_from_fun, data=adjsurv_0,
                     est=est, interpolation=interpolation,
                     FUN.VALUE=numeric(1))
    surv_1 <- vapply(times, read_from_fun, data=adjsurv_1,
                     est=est, interpolation=interpolation,
                     FUN.VALUE=numeric(1))

    if (conf_int) {
      se_0 <- vapply(times, read_from_step_function, data=adjsurv_0,
                     est="se", FUN.VALUE=numeric(1))
      se_1 <- vapply(times, read_from_step_function, data=adjsurv_1,
                     est="se", FUN.VALUE=numeric(1))
    }
  }

  surv_diff <- surv_0 - surv_1

  diff_dat <- data.frame(time=times)
  diff_dat[, est] <- surv_diff

  if (conf_int) {
    # calculate confidence interval from pooled SE
    diff_se <- sqrt(se_0^2 + se_1^2)
    diff_ci <- confint_surv(surv=surv_diff, se=diff_se, conf_level=conf_level,
                            conf_type="plain")
    # put together
    diff_dat$se <- diff_se
    diff_dat$ci_lower <- diff_ci$left
    diff_dat$ci_upper <- diff_ci$right
  }
  return(diff_dat)
}
