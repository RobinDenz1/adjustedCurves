####################### Helper functions for R-Package #########################

## estimate iptw weights
get_iptw_weights <- function(data, treatment_model, weight_method, ...) {

  if (inherits(treatment_model, "formula")) {

    args <- list(formula=treatment_model, data=data,
                 method=weight_method,
                 estimand="ATE")
    weights <- do.call(WeightIt::weightit, c(args, ...))$weights

    # TODO: allow all objects with predict method
  } else if (inherits(treatment_model, "glm")) {

    ps <- predict(treatment_model, newdata=data, type="response")
    weights <- ifelse(data[, variable]==1, 1/ps, 1/(1-ps))

  }
  return(weights)
}

## relevant event times
get_times <- function(data, variable, ev_time, event, type,
                      custom=NULL) {

  levels <- sort(unique(data[, variable]))

  if (type=="events_overall") {

    time <- sort(unique(data[, ev_time][data[, event]==1]))
    times <- list()
    for (i in 1:length(levels)) {
      times[[toString(levels[i])]] <- time
    }
    times[["overall"]] <- time

  } else if (type=="events_group") {

    times <- list()
    for (i in 1:length(levels)) {
      time <- sort(unique(data[, ev_time][data[, event]==1 &
                                            data[, variable]==levels[i]]))
      times[[toString(levels[i])]] <- time
    }
    times[["overall"]] <- sort(unique(data[, ev_time][data[, event]==1]))

  } else if (type=="vector") {

    times <- list()
    for (i in 1:length(levels)) {
      times[[toString(levels[i])]] <- custom
    }
    times[["overall"]] <- custom

  } else if (type=="list") {
    # some input checks i guess
  }

  return(times)

}

# calculate CI from sd
confint_surv <- function(surv, sd, n, alpha) {
  error <- qnorm(1-(alpha/2)) * sd / sqrt(n)
  left <- surv - error
  right <- surv + error
  return(list(left=left, right=right))
}

## throw errors when inputs don't make sense
# TODO: - doesn't care about defaults in functions sadly
#       - should also check if needed arguments are supplied
check_inputs_adjusted_surv <- function(data, variable, ev_time, event, method,
                                       sd, times, ...) {
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
  }

  # Here: check if times input is correct if method requires it
  if (method %in% c("iptw_pseudo", "direct", "aiptw", "direct_pseudo",
                    "aiptw_pseudo", "el", "tmle", "ostmle")) {
    # TODO: time checks
  } else {
    if (!is.null(obj$times)) {
      warning("Object 'times' is not defined for method='",
              method, "' and will be ignored.")
    }
  }

  # Direct Pseudo, AIPTW Pseudo
  if (method=="direct_pseudo" | method=="aiptw_pseudo" | method=="iptw_pseudo") {
    require("geese")
    require("prodlim")

    if (!is.character(obj$outcome_vars)) {
      stop("'outcome_vars' should be a character vector of column names",
           "in 'data', used to model the outcome mechanism.")
    } else if (!obj$type_time %in% c("factor", "bs", "ns")) {
      stop("'type_time' should be either 'factor', 'bs' or 'ns'.")
    } else if (!is.logical(obj$na.rm)) {
      stop("Argument 'na.rm' must be either TRUE or FALSE.")
    }

    if (method=="aiptw_pseudo" | method=="iptw_pseudo") {
      if (!(is.numeric(obj$treatment_model) | inherits(obj$treatment_model, "glm"))) {
        stop("Argument 'treatment_model' must be either a glm object or a",
             " numeric vector of propensity scores.")
      }
    }
  # TMLE
  } else if (method=="tmle") {
    require("survtmle")

    if (!all(obj$times$overall==floor(obj$times$overall))) {
      stop("Only integer time is allowed when using method='tmle'.")
    }
  # OSTMLE
  } else if (method=="ostmle") {
    require("MOSS")

    if (!all(obj$times$overall==floor(obj$times$overall))) {
      stop("Only integer time is allowed when using method='ostmle'.")
    }
  # Empirical Likelihood
  } else if (method=="el") {
    require("adjKMtest")

    if (!is.character(obj$treatment_vars)) {
      stop("'treatment_vars' should be a character vector of column names",
           "in 'data', used to model the outcome mechanism.")
    } else if (!all(obj$treatment_vars %in% colnames(data))) {
      stop("'treatment_vars' should be a character vector of column names",
           "in 'data', used to model the outcome mechanism.")
    } else if (!obj$moment %in% c("first", "second")) {
      stop("Argument 'moment' must be either 'first' or 'second'.")
    } else if (!is.logical(obj$standardize)) {
      stop("Argument 'standardize' must be either TRUE or FALSE.")
    } # TODO: check for dichotomous variables, warn when 0, 1
  }
}
