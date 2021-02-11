## for adjustedsurv function
check_inputs_adjustedsurv <- function(data, variable, ev_time, event, method,
                                      conf_int, conf_level, times, bootstrap,
                                      n_boot, na.rm, ...) {
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
  # conf_int
  } else if (!is.logical(conf_int)) {
    stop("'conf_int' must be either TRUE or FALSE.")
  } else if (!is.logical(na.rm)) {
    stop("'na.rm' must be either TRUE or FALSE.")
  } else if (!is.numeric(conf_level)) {
    stop("'conf_level' must be a number < 1 and > 0.")
  } else if (conf_level >= 1 | conf_level <= 0) {
    stop("'conf_level' must be a number < 1 and > 0.")
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

    if (method=="aiptw") {
      if ((!"censoring_model" %in% names(obj)) &
          (!"treatment_model" %in% names(obj)) &
          (!"outcome_model" %in% names(obj))) {
        stop("At least one of 'treatment_model', 'outcome_model' and ",
             "'censoring_model' needs to be specified, see details.")
      }
    }

    # TMLE
  } else if (method=="tmle") {
    requireNamespace("survtmle")

    if (!is.null(obj$times)) {
      if (!all(obj$times==floor(obj$times))) {
        stop("Only integer time is allowed when using method='tmle'.")
      }
    }
    # OSTMLE
  } else if (method=="ostmle") {
    requireNamespace("MOSS")

    if (!is.null(obj$times)) {
      if (!all(obj$times==floor(obj$times))) {
        stop("Only integer time is allowed when using method='ostmle'.")
      }
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
    } else if ("moment" %in% names(obj)) {

      if (!obj$moment %in% c("first", "second")) {
        stop("Argument 'moment' must be either 'first' or 'second'.")
      }

    } else if ("standardize" %in% names(obj)) {
      if (!is.logical(obj$standardize)) {
        stop("Argument 'standardize' must be either TRUE or FALSE.")
      }
    } else if (conf_int) {
      warning("No asymptotic or exact variance approximations exist for ",
              "method='el'. Use bootstrap=TRUE to get bootstrap estimates.")
    }
     # TODO: check for dichotomous variables, warn when 0, 1
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

## for sim_confounded_surv function
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

## for adjusted_rmst function
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

## for adjustedsurv_test
check_inputs_adj_test <- function(adjsurv, from, to) {

  if (!inherits(adjsurv, "adjustedsurv")) {
    stop("'adjsurv' must be an 'adjustedsurv' object, created using the ",
         "adjustedsurv function.")
  } else if (is.null(adjsurv$boot_data)) {
    stop("Can only perform a significance test if bootstrapping was ",
         "performed (bootstrap=TRUE in adjsutedsurv() call).")
  } else if (adjsurv$categorical) {
    stop("This function currently only supports a test of two survival curves.")
  }

}
