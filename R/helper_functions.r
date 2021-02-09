####################### Helper functions for R-Package #########################

## estimate iptw weights
get_iptw_weights <- function(data, treatment_model, weight_method,
                             variable, ...) {

  if (inherits(treatment_model, "formula")) {

    args <- list(formula=treatment_model, data=data,
                 method=weight_method,
                 estimand="ATE")
    weights <- do.call(WeightIt::weightit, c(args, ...))$weights

  } else if (inherits(treatment_model, "glm")) {

    ps <- stats::predict.glm(treatment_model, newdata=data, type="response")
    weights <- ifelse(data[, variable]==1, 1/ps, 1/(1-ps))

  } else if (inherits(treatment_model, "multinom")) {

    predict.multinom <- utils::getFromNamespace("predict.multinom", "nnet")
    ps <- predict.multinom(treatment_model, newdata=data, type="probs")

    weights <- rep(0, nrow(data))
    for (i in levels(data[,variable])) {
      weights[data[,variable] == i] <- 1/ps[data[,variable] == i, i]
    }

  }
  return(weights)
}

## calculate CI from sd
# TODO: This can't be correct ...
confint_surv <- function(surv, sd, n, alpha, conf_type="plain") {
  if (conf_type=="plain") {
    error <- stats::qnorm(1-(alpha/2)) * sd / sqrt(n)
    left <- surv - error
    right <- surv + error
  } else if (conf_type=="log") {
    error <- stats::qnorm(1-(alpha/2)) * sd / sqrt(n)
    left <- surv * exp(-error)
    right <- surv * exp(error)
  }

  return(list(left=left, right=right))
}

## simulate survival time according to Bender et al. (2005)
sim_surv_time <- function(row, betas, dist, lambda,
                          gamma) {
  U <- stats::runif(1, min=0, max=1)
  eff <- sum(row * betas)

  if (dist=="weibull") {
    surv_time <- (-(log(U)/(lambda*exp(eff))))^(1/gamma)
  } else if (dist=="exponential") {
    surv_time <- -(log(U)/(lambda*exp(eff)))
  }
  return(surv_time)
}

## takes a value x at which to read from the step function
## and step function data from which to read it
read_from_step_function <- function(x, step_data) {
  check <- step_data[which(step_data$time <= x),]
  if (nrow(check)==0) {
    val <- 1
  } else {
    val <- check$surv[which(check$time==max(check$time))][1]
  }
  return(val)
}

## calculate difference between two step functions
## according to some transformation function
exact_stepfun_difference <- function(adjsurv, times, max_t) {

  times <- times[times<=max_t]

  levs <- unique(adjsurv$group)
  adjsurv_0 <- adjsurv[which(adjsurv$group==levs[1]),]
  adjsurv_1 <- adjsurv[which(adjsurv$group==levs[2]),]

  surv_0 <- sapply(times, read_from_step_function, step_data=adjsurv_0)
  surv_1 <- sapply(times, read_from_step_function, step_data=adjsurv_1)
  surv_diff <- surv_1 - surv_0

  diff_dat <- data.frame(time=times, surv=surv_diff)

  return(diff_dat)
}

## calculate exact integral under step function
# 'stepfun' needs to be a data.frame with columns 'time' and 'surv',
# sorted by time with no duplicates in time
exact_stepfun_integral <- function(stepfun) {
  integral <- 0
  for (i in 1:(length(stepfun$time)-1)) {
    x1 <- stepfun$time[i]
    x2 <- stepfun$time[i+1]
    y <- stepfun$surv[i]
    rect_area <- (x2 - x1) * y
    integral <- integral + rect_area
  }
  return(integral)
}

## throw errors when inputs don't make sense
check_inputs_adjustedsurv <- function(data, variable, ev_time, event, method,
                                      sd, times, bootstrap, n_boot, na.rm, ...) {
  requireNamespace("survival")

  obj <- list(...)

  if (!inherits(data, "data.frame")) {
    stop("'data' argument must be a data.frame object.")
  # needed variables
  } else if (!is.character(variable) | !is.character(ev_time) |
             !is.character(event) | !is.character(method)) {
    stop("Arguments 'variable', 'ev_time', 'event' and 'method' must be ",
         "character strings, specifying variables in 'data'.")
  } else if (!variable %in% colnames(data)) {
    stop(variable, " is not a valid column name in 'data'.")
  } else if (!ev_time %in% colnames(data)) {
    stop(ev_time, " is not a valid column name in 'data'.")
  } else if (!event %in% colnames(data)) {
    stop(event, " is not a valid column name in 'data'.")
  # method
  } else if (!method %in% c("km", "iptw_km", "iptw_cox", "iptw_pseudo",
                            "direct", "direct_pseudo", "aiptw_pseudo",
                            "aiptw", "tmle", "ostmle", "matching", "el")) {
    stop("Method '", method, "' is undefined. See documentation for ",
         "details on available methods.")
  # sd
  } else if (!is.logical(sd)) {
    stop("'sd' must be either TRUE or FALSE.")
  } else if (!is.logical(na.rm)) {
    stop("'na.rm' must be either TRUE or FALSE.")
  }

  # Check if the group variable has the right format
  if (method %in% c("matching", "el", "tmle", "ostmle") &
      is.factor(data[,variable])) {
    stop("The column in 'data' specified by 'variable' needs to be ",
         "a dichotomous integer variable if method='", method, "'.")
  }

  if (!method %in% c("matching", "el", "tmle", "ostmle") &
      !is.factor(data[,variable])) {
    stop("The column in 'data' specified by 'variable' needs to be ",
         "a factor variable if method='", method, "'.")
  }

  # Check if categorical should be allowed
  if (length(unique(data[,variable])) > 2 &
      method %in% c("matching", "el", "tmle", "ostmle", "aiptw")) {
    stop("Categorical treatments are currently not supported for ",
         "method='", method, "'.")
  }

  # Here: check if times input is correct if method requires it
  if (method %in% c("iptw_pseudo", "direct", "aiptw", "direct_pseudo",
                    "aiptw_pseudo", "el", "tmle", "ostmle")) {
    if (!is.numeric(times) & !is.null(times)) {
      stop("'times' must be a numeric vector.")
    }
  } else {
    if (!is.null(obj$times)) {
      warning("Object 'times' is not defined for method='",
              method, "' and will be ignored.")
    }
  }

  # Direct Pseudo, AIPTW Pseudo
  if (method=="direct_pseudo" | method=="aiptw_pseudo" | method=="iptw_pseudo") {
    requireNamespace("geepack")
    requireNamespace("prodlim")

    if ("outcome_vars" %in% names(obj)) {
      if (!is.character(obj$outcome_vars)) {
        stop("'outcome_vars' should be a character vector of column names ",
             "in 'data', used to model the outcome mechanism.")
      }
    }

    if ("type_time" %in% names(obj)) {
      if (!obj$type_time %in% c("factor", "bs", "ns")) {
        stop("'type_time' should be either 'factor', 'bs' or 'ns'.")
      }
    }

    if ("treatment_model" %in% names(obj)) {
      if (method=="aiptw_pseudo" | method=="iptw_pseudo") {
        if (!(is.numeric(obj$treatment_model) |
              inherits(obj$treatment_model, "glm") |
              inherits(obj$treatment_model, "multinom"))) {
          stop("Argument 'treatment_model' must be either a glm object or a",
               " numeric vector of propensity scores.")
        }
      }
    } else {
      stop("Argument 'treatment_model' is missing with no standard value.")
    }
  # TMLE
  } else if (method=="tmle") {
    requireNamespace("survtmle")

    if (!all(obj$times==floor(obj$times))) {
      stop("Only integer time is allowed when using method='tmle'.")
    }
  # OSTMLE
  } else if (method=="ostmle") {
    requireNamespace("MOSS")

    if (!all(obj$times==floor(obj$times))) {
        stop("Only integer time is allowed when using method='ostmle'.")
    }
  # Empirical Likelihood
  } else if (method=="el") {
    requireNamespace("adjKMtest")

    if (!is.character(obj$treatment_vars)) {
      stop("'treatment_vars' should be a character vector of column names ",
           "in 'data', used to model the outcome mechanism.")
    } else if (!all(obj$treatment_vars %in% colnames(data))) {
      stop("'treatment_vars' should be a character vector of column names ",
           "in 'data', used to model the outcome mechanism.")
    } else if (!obj$moment %in% c("first", "second")) {
      stop("Argument 'moment' must be either 'first' or 'second'.")
    } else if (!is.logical(obj$standardize)) {
      stop("Argument 'standardize' must be either TRUE or FALSE.")
    } # TODO: check for dichotomous variables, warn when 0, 1
  }

  # bootstrapping
  if (bootstrap) {

    if (is.numeric(obj$treatment_model)) {
      stop("'treatment_model' needs to be an actual model that can be ",
           "refit when using bootstrap=TRUE.")
    } else if (is.numeric(obj$outcome_model)) {
      stop("'outcome_model' needs to be an actual model that can be ",
           "refit when using bootstrap=TRUE.")
    }

  }
}

## throw error when inputs don't make sense
check_inputs_sim_fun <- function(n, lcovars, outcome_betas, surv_dist,
                                 gamma, lambda, treatment_betas,
                                 intercept, gtol, cens_fun, cens_args,
                                 max_t, group_beta) {

  if (!is.numeric(n)) {
    stop("'n' must be a positive integer.")
  } else if (floor(n)!=n) {
    stop("'n' must be a positive integer, not a double.")
  } else if(!is.character(surv_dist)) {
    stop("'surv_dist' must be either 'weibull' or 'exponential'.")
  } else if (!surv_dist %in% c("weibull", "exponential")) {
    stop("'surv_dist' must be either 'weibull' or 'exponential'.")
  } else if (!is.numeric(gamma)) {
    stop("'gamma' must be a number. See details.")
  } else if (!is.numeric(lambda)) {
    stop("'lambda' must be a number. See details.")
  } else if (!is.numeric(intercept)) {
    stop("'intercept' must be a number. See details.")
  } else if (!is.numeric(gtol)) {
    stop("'gtol' must be a number. See details.")
  } else if (gtol > 1 | gtol < 0) {
    stop("'gtol' must be <= 1 and >= 0. See details.")
  } else if (!is.function(cens_fun)) {
    stop("'cens_fun' must be a function with the argument 'n'.")
  } else if (!is.list(cens_args)) {
    stop("'cens_args' must be a named list of arguments to be passed to",
         " 'cens_fun'.")
  } else if (!is.numeric(max_t)) {
    stop("'max_t' must be a number.")
  } else if (max_t <= 0) {
    stop("'max_t' must be bigger than zero.")
  } else if (!is.numeric(group_beta)) {
    stop("'group_beta' must be a number.")
  }

  if (!((is.null(lcovars) & is.null(outcome_betas) & is.null(treatment_betas)) |
      (is.list(lcovars) & is.numeric(outcome_betas) & is.numeric(treatment_betas)))) {
    stop("'lcovars', 'outcome_betas' and 'treatment_betas' can either all ",
         "be NULL to use default values, or have to be specified.")
  }

  if (!is.null(lcovars)) {
    if (is.null(names(lcovars))) {
      stop("Elements in the 'lcovars' list must be named.")
    }
  }

  if (!is.null(outcome_betas)) {
    if (is.null(names(outcome_betas))) {
      stop("Elements in the 'outcome_betas' vector must be named.")
    }
  }

  if (!is.null(treatment_betas)) {
    if (is.null(names(treatment_betas))) {
      stop("Elements in the 'treatment_betas' vector must be named.")
    }
  }

  if (!is.null(lcovars)) {
    if (!(names(lcovars) == names(outcome_betas) &
        names(treatment_betas) == names(lcovars) &
        names(outcome_betas) == names(treatment_betas))) {
      stop("The names of the objects in 'lcovars', 'outcome_betas' and ",
           " 'treatment_betas' must be the same.")
    }
  }

  if (!is.null(lcovars)) {
    for (i in 1:length(lcovars)) {
      if (!lcovars[[i]][1] %in% c("rbinom", "rnorm", "runif")) {
        stop("The first element of every vector in 'lcovars' must be either ",
             "'rbinom', 'rnorm' or 'runif', not ", lcovars[[i]][1])
      }
    }
  }
}

## throw error when inputs don't make sense
check_inputs_adj_rmst <- function(adjsurv, from, to, use_boot) {

  if (!is.numeric(from) | !is.numeric(to)) {
    stop("'from' and 'to' must be numbers.")
  } else if (!inherits(adjsurv, "adjustedsurv")) {
    stop("'adjsurv' must be an 'adjustedsurv' object created using ",
         "the 'adjustedsurv()' function.")
  } else if (from >= to) {
    stop("'from' must be smaller than 'to'.")
  } else if (use_boot & is.null(adjsurv$boot_data)) {
    warning("Cannot use bootstrapped estimates because they were not estimated.",
            " Need 'bootstrap=TRUE' in 'adjustedsurv' function call.")
  }

}
