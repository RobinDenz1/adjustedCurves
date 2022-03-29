
## takes a value x at which to read from the step function
## and step function data from which to read it
read_from_step_function <- function(x, data, est="surv") {

  # keep only data with non-missing est
  data <- data[which(!is.na(data[, est])), ]

  # no extrapolation
  if (x > max(data$time)) {
    return(NA)
  }

  # otherwise get value
  check <- data[which(data$time <= x), ]
  if (nrow(check)==0) {
    if (est=="surv") {
      val <- 1
    } else if (est=="cif") {
      val <- 0
    } else {
      val <- NA
    }
  } else {
    val <- check[, est][which(check$time==max(check$time))][1]
  }
  return(val)
}

## takes a value x at which to read from the linear function
## and step function data from which to read it
read_from_linear_function <- function(x, data, est="surv") {

  time_vec <- data$time
  est_vec <- data[, est]

  # add zero to beginning if necessary
  if (est=="surv") {
    val <- suppressWarnings(stats::approx(x=time_vec, y=est_vec, xout=x,
                                          yleft=0)$y)
  } else if (est=="cif") {
    val <- suppressWarnings(stats::approx(x=time_vec, y=est_vec, xout=x,
                                          yleft=1)$y)
  } else {
    val <- suppressWarnings(stats::approx(x=time_vec, y=est_vec, xout=x)$y)
  }
  return(val)
}

## calculate integral of a function using the trapezoid rule
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

## calculate exact integral under step function
# 'stepfun' needs to be a data.frame (!) with columns 'time' and 'est',
# sorted by time with no duplicates in time
exact_stepfun_integral <- function(stepfun, from, to, est="surv") {

  # constrain step function end
  latest <- read_from_step_function(to, data=stepfun, est=est)
  stepfun <- stepfun[stepfun$time <= to, ]

  if (!to %in% stepfun$time) {
    temp <- data.frame(time=to)
    temp[, est] <- latest
    stepfun <- rbind(stepfun, temp)
  }

  # constrain step function beginning
  if (from != 0) {
    earliest <- read_from_step_function(from, data=stepfun, est=est)
    stepfun <- stepfun[stepfun$time >= from, ]

    if (!from %in% stepfun$time) {
      temp <- data.frame(time=from)
      temp[, est] <- earliest
      stepfun <- rbind(temp, stepfun)
    }
  }

  # when there are unknown survival times, return NA
  # only necessary when the last time is NA, since my
  # algorithm technically still works in that case, but makes no sense
  if (anyNA(stepfun)) {
    return(NA)
  }

  # calculate exact integral
  integral <- 0
  for (i in seq_len((length(stepfun$time)-1))) {
    x1 <- stepfun$time[i]
    x2 <- stepfun$time[i+1]
    y <- stepfun[, est][i]
    rect_area <- (x2 - x1) * y
    integral <- integral + rect_area
  }
  return(integral)
}

## calculate difference between two functions at arbitrary time values
## using either linear or step-function interpolation
difference_function <- function(adj, times, est="surv", type="steps") {

  levs <- levels(adj$group)
  adjsurv_0 <- adj[which(adj$group==levs[1]), ]
  adjsurv_1 <- adj[which(adj$group==levs[2]), ]

  if (type=="steps") {
    read_fun <- read_from_step_function
  } else if (type=="linear") {
    read_fun <- read_from_linear_function
  }

  if (nrow(adjsurv_0)==nrow(adjsurv_1) && all(adjsurv_0$time==adjsurv_1$time) &&
      type=="steps") {
    surv_0 <- adjsurv_0[, est]
    surv_1 <- adjsurv_1[, est]

  } else {
    surv_0 <- vapply(times, read_fun, data=adjsurv_0,
                     est=est, FUN.VALUE=numeric(1))
    surv_1 <- vapply(times, read_fun, data=adjsurv_1,
                     est=est, FUN.VALUE=numeric(1))
  }

  surv_diff <- surv_0 - surv_1

  diff_dat <- data.frame(time=times)
  diff_dat[, est] <- surv_diff

  return(diff_dat)
}

## function to calculate the integral of the difference of two functions,
## using either linear or step-function interpolation
difference_integral <- function(adj, from, to, type="steps", est="surv",
                                subdivisions=1000) {
  if (type=="linear") {
    times <- seq(from, to, (to-from)/subdivisions)
    diff_dat <- difference_function(adj=adj, times=times, est=est,
                                    type="linear")
    area <- trapezoid_integral(x=diff_dat$time, y=diff_dat[, est])
  } else if (type=="steps") {
    times <- sort(unique(adj$time))
    diff_dat <- difference_function(adj=adj, times=times, est=est,
                                    type="steps")
    area <- exact_stepfun_integral(stepfun=diff_dat, from=from, to=to,
                                   est=est)
  }

  output <- list(times=times,
                 diff_dat=diff_dat,
                 area=area)
  return(output)
}
