## for adjustedsurv function
check_inputs_adjustedsurv <- function(data, variable, ev_time, event, method,
                                      conf_int, conf_level, times, bootstrap,
                                      n_boot, na.action, ...) {
  requireNamespace("survival")

  obj <- list(...)

  if (!inherits(data, c("data.frame", "mids"))) {
    stop("'data' argument must be a data.frame or mids object.")
  # needed variables
  } else if (!is.character(variable) | !is.character(ev_time) |
             !is.character(event) | !is.character(method)) {
    stop("Arguments 'variable', 'ev_time', 'event' and 'method' must be ",
         "character strings, specifying variables in 'data'.")
  } else if (!method %in% c("km", "iptw_km", "iptw_cox", "iptw_pseudo",
                            "direct", "direct_pseudo", "aiptw_pseudo",
                            "aiptw", "tmle", "ostmle", "matching", "emp_lik")) {
    stop("Method '", method, "' is undefined. See documentation for ",
         "details on available methods.")
  # conf_int
  } else if (!is.logical(conf_int)) {
    stop("'conf_int' must be either TRUE or FALSE.")
  } else if (!is.numeric(conf_level)) {
    stop("'conf_level' must be a number < 1 and > 0.")
  } else if (conf_level >= 1 | conf_level <= 0) {
    stop("'conf_level' must be a number < 1 and > 0.")
  }

  if (inherits(data, "data.frame")) {
    if (!variable %in% colnames(data)) {
      stop(variable, " is not a valid column name in 'data'.")
    } else if (!ev_time %in% colnames(data)) {
      stop(ev_time, " is not a valid column name in 'data'.")
    } else if (!event %in% colnames(data)) {
      stop(event, " is not a valid column name in 'data'.")
    }

    # Check if categorical should be allowed
    if (length(unique(data[,variable])) > 2 &
        method %in% c("matching", "emp_lik", "tmle", "ostmle",
                      "aiptw")) {
      stop("Categorical treatments are currently not supported for ",
           "method='", method, "'.")
    }

    # Check if the group variable has the right format
    if (method %in% c("matching", "emp_lik", "tmle", "ostmle") &
        !is.factor(data[,variable]) & !is.numeric(data[,variable])) {
      stop("The column in 'data' specified by 'variable' needs to be ",
           "a factor or a dichotomous integer variable if method='", method, "'.")
    }

    if (!method %in% c("matching", "emp_lik", "tmle", "ostmle") &
        !is.factor(data[,variable])) {
      stop("The column in 'data' specified by 'variable' needs to be ",
           "a factor variable if method='", method, "'.")
    }

    # No extrapolation
    if (!is.null(times)) {
      if (max(times) > max(data[,ev_time])) {
        stop("Values in '", ev_time, "' must be smaller than max(data[,ev_time]).",
             " No extrapolation allowed.")
      }
    }
  } else {
    if (!is.null(obj$treatment_model) & !inherits(obj$treatment_model, c("mira", "formula"))) {
      stop("When using multiple imputation, mira objects need to be supplied",
           " to 'treatment_model' instead of single models. See documentation.")
    }
    if (!is.null(obj$outcome_model) & !inherits(obj$outcome_model, "mira")) {
      stop("When using multiple imputation, mira objects need to be supplied",
           " to 'outcome_model' instead of single models. See documentation.")
    }
    if (!is.null(obj$censoring_model) & !inherits(obj$censoring_model, "mira")) {
      stop("When using multiple imputation, mira objects need to be supplied",
           " to 'censoring_model' instead of single models. See documentation.")
    }
  }

  # Here: check if times input is correct if method requires it
  if (method %in% c("iptw_pseudo", "direct", "aiptw", "direct_pseudo",
                    "aiptw_pseudo", "emp_lik", "tmle", "ostmle")) {
    if (!is.numeric(times) & !is.null(times)) {
      stop("'times' must be a numeric vector.")
    }
  } else {
    if (!is.numeric(times) & !is.null(times)) {
      warning("Object 'times' must be a numeric vector or NULL if method='",
              method, "'.")
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
      if (method=="aiptw_pseudo") {
        if (!(is.numeric(obj$treatment_model) |
              inherits(obj$treatment_model, c("glm", "multinom", "mira")))) {
          stop("Argument 'treatment_model' must be one of: 'glm',",
          " 'multinom', 'mira' or a numeric vector of propensity scores.")
        }
      }
    }

    if (method %in% c("aiptw_pseudo", "direct_pseudo") & !is.null(times)) {
      if (length(times)==1) {
        stop("'geese' models require at least two distinct time points. ",
             "Add more points in time to 'times' and try again.")
      }
      if (!is.null(obj$spline_df)) {
        if (obj$spline_df > length(times)) {
          warning("'spline_df' > len(times) might lead to problems.")
        }
      } else if (!is.null(obj$type_time)) {
        if (5 > length(times) & obj$type_time!="factor") {
          warning("'spline_df' > length(times)=5 might lead to problems when",
                  " type_time!='factor'.")
        }
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

    if (!is.null(obj$times)) {
      if (!all(obj$times==floor(obj$times))) {
        stop("Only integer time is allowed when using method='ostmle'.")
      }
    }
  # Empirical Likelihood
  } else if (method=="emp_lik") {
    requireNamespace("MASS")

    if (!is.character(obj$treatment_vars)) {
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
    } else if (inherits(data, "data.frame")) {

      if (!all(obj$treatment_vars %in% colnames(data))) {
        stop("'treatment_vars' should be a character vector of column names ",
             "in 'data', used to model the outcome mechanism.")
      }

      for (col in obj$treatment_vars) {
        if (paste0(unique(data[,col]), collapse="") %in% c("10", "01")) {
          warning("Dichotomous variables coded with 0 and 1 found in 'treatment_vars'.",
                  " Consider recoding to -1 and 1 to avoid estimation problems.")
        }
      }

    }
  # Matching
  } else if (method=="matching") {
    requireNamespace("Matching")

    if (bootstrap) {
      warning("Bootstrapping generally doesn't produce unbiased variance",
              " estimates with matching estimators. Use with caution. ",
              "See ?surv_matching.")
    } else if (is.numeric(obj$treatment_model)) {
      if (any(obj$treatment_model > 1) | any(obj$treatment_model < 0)) {
        stop("Propensity Scores > 1 or < 0 not allowed. Perhaps you supplied ",
             "weights on accident?")
      }
    }
  # IPTW KM
  } else if (method=="iptw_km") {
    if (conf_int & inherits(obj$treatment_model, "formula")) {
      warning("Approximate confidence intervals currently not supported in ",
              "method='iptw_km' when 'treatment_model' is not a 'glm' or 'multinom'",
              " object.")
    }
  # AIPTW
  } else if (method=="aiptw") {
    if ((!"censoring_model" %in% names(obj)) &
        (!"treatment_model" %in% names(obj)) &
        (!"outcome_model" %in% names(obj))) {
      stop("At least one of 'treatment_model', 'outcome_model' and ",
           "'censoring_model' needs to be specified, see details.")
    }
  } else if (method=="direct") {
    if (!"outcome_model" %in% names(obj)) {
      stop("Argument 'outcome_model' must be specified when using method='direct'.")
    }
  }

  # bootstrapping
  if (bootstrap) {

    if (is.numeric(obj$treatment_model)) {
      stop("'treatment_model' needs to be a model that can be ",
           "refit or a formula object when using bootstrap=TRUE.")
    } else if (is.numeric(obj$outcome_model)) {
      stop("'outcome_model' needs to be an actual model that can be ",
           "refit when using bootstrap=TRUE.")
    }

  }

  # asymptotic variance calculations
  if (conf_int) {
    if (method %in% c("emp_lik", "matching")) {
      warning("Asymptotic or exact variance calculations are currently",
              " not available for method='", method, "'. Use bootstrap=TRUE",
              "to get bootstrap estimates.")
    } else if (method=="direct_pseudo" & !is.null(obj$model_type)) {
      if (obj$model_type != "lm") {
        warning("Asymptotic variance calculations for method='direct_pseudo' ",
                "can only be calculated with model_type='lm'.")
      }
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
  } else if (!is.function(cens_fun) & !is.null(cens_fun)) {
    stop("'cens_fun' must be a function with the argument 'n' or NULL.")
  } else if (!is.list(cens_args) & !is.null(cens_args)) {
    stop("'cens_args' must be a named list of arguments to be passed to",
         " 'cens_fun' or NULL.")
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
    if (!(all(names(lcovars) == names(outcome_betas)) &
          all(names(treatment_betas) == names(lcovars)) &
          all(names(outcome_betas) == names(treatment_betas)))) {
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

  if (!(inherits(adjsurv, "adjustedsurv") | inherits(adjsurv, "adjustedcif"))) {
    stop("'adjsurv' must be an 'adjustedsurv' object, created using the ",
         "adjustedsurv function.")
  } else if (is.null(adjsurv$boot_data)) {
    stop("Can only perform a significance test if bootstrapping was ",
         "performed (bootstrap=TRUE in adjustedsurv/adjustedcif call).")
  } else if (!is.numeric(from)) {
    stop("'from' must be a number >= 0.")
  } else if (from < 0) {
    stop("'from' must be a number >= 0.")
  } else if (!is.numeric(to)) {
    stop("'to' must be a number >= 0.")
  } else if (to <= from) {
    stop("'to' must be greater than 'from'.")
  }

  if (class(adjsurv)=="adjustedsurv") {
    if (to > max(adjsurv$adjsurv$time)) {
      stop("'to' can not be greater than the latest observed time.")
    }
  } else {
    if (to > max(adjsurv$adjcif$time)) {
      stop("'to' can not be greater than the latest observed time.")
    }
  }

}

## for adjustedcif
check_inputs_adjustedcif <- function(data, variable, ev_time, event, method,
                                     conf_int, conf_level, times, bootstrap,
                                     n_boot, cause=cause, na.action, ...) {
  requireNamespace("survival")

  obj <- list(...)

  if (!inherits(data, c("data.frame", "mids"))) {
    stop("'data' argument must be a data.frame or mids object.")
  # needed variables
  } else if (!is.character(variable) | !is.character(ev_time) |
             !is.character(event) | !is.character(method)) {
    stop("Arguments 'variable', 'ev_time', 'event' and 'method' must be ",
         "character strings, specifying variables in 'data'.")
  }

  if (inherits(data, "data.frame")) {
    # Check input variables format
    if (!variable %in% colnames(data)) {
      stop(variable, " is not a valid column name in 'data'.")
    } else if (!ev_time %in% colnames(data)) {
      stop(ev_time, " is not a valid column name in 'data'.")
    } else if (!event %in% colnames(data)) {
      stop(event, " is not a valid column name in 'data'.")
    }

    # Check if the group variable has the right format
    if (method %in% c("matching", "tmle") &
        !is.factor(data[,variable]) & !is.numeric(data[,variable])) {
      stop("The column in 'data' specified by 'variable' needs to be ",
           "a dichotomous integer variable or a factor variable if method='",
           method, "'.")
    }

    if (!method %in% c("matching", "tmle") &
        !is.factor(data[,variable])) {
      stop("The column in 'data' specified by 'variable' needs to be ",
           "a factor variable if method='", method, "'.")
    }

    # Check if categorical should be allowed
    if (length(unique(data[,variable])) > 2 &
        method %in% c("matching", "tmle", "aiptw")) {
      stop("Categorical treatments are currently not supported for ",
           "method='", method, "'.")
    }

  }

  if (!method %in% c("aalen_johansen", "iptw", "iptw_pseudo", "direct",
                            "direct_pseudo", "aiptw_pseudo",
                            "aiptw", "tmle", "matching")) {
    stop("Method '", method, "' is undefined. See documentation for ",
         "details on available methods.")
  # conf_int
  } else if (!is.logical(conf_int)) {
    stop("'conf_int' must be either TRUE or FALSE.")
  } else if (!is.numeric(conf_level)) {
    stop("'conf_level' must be a number < 1 and > 0.")
  } else if (conf_level >= 1 | conf_level <= 0) {
    stop("'conf_level' must be a number < 1 and > 0.")
  # cause
  } else if (!is.numeric(cause)) {
    stop("'cause' must be a number specifying the cause of interest in ",
         "the column specified with 'event'.")
  # time
  } else if (!is.numeric(times) & !is.null(times)) {
    stop("'times' must be a numeric vector or NULL.")
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
      if (method=="aiptw_pseudo") {
        if (!(is.numeric(obj$treatment_model) |
              inherits(obj$treatment_model, c("glm", "multinom", "mira")))) {
          stop("Argument 'treatment_model' must be one of glm, multinom or mira",
               " or a numeric vector of propensity scores.")
        }
      }
    }

    if (method %in% c("aiptw_pseudo", "direct_pseudo") & !is.null(times)) {
      if (length(times)==1) {
        stop("'geese' models require at least two distinct time points. ",
             "Add more points in time to 'times' and run again.")
      }
      if (!is.null(obj$spline_df)) {
        if (obj$spline_df > length(times)) {
          warning("'spline_df' > len(times) might lead to problems.")
        }
      } else if (!is.null(obj$type_time)) {
        if (5 > length(times) & obj$type_time!="factor") {
          warning("'spline_df' > len(times) might lead to problems when",
                  " type_time!='factor'.")
        }
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
  # Matching
  } else if (method=="matching") {
    requireNamespace("Matching")

    if (bootstrap) {
      warning("Bootstrapping generally doesn't produce unbiased variance",
              " estimates with matching estimators. Use with caution. ",
              "See ?surv_method_matching.")
    } else if (is.numeric(obj$treatment_model)) {
      if (any(obj$treatment_model > 1) | any(obj$treatment_model < 0)) {
        stop("Propensity Scores > 1 or < 0 not allowed. Perhaps you supplied ",
             "weights on accident?")
      }
    }

  } else if (method=="aalen_johansen") {
    requireNamespace("cmprsk")
  }

  # bootstrapping
  if (bootstrap) {

    if (is.numeric(obj$treatment_model)) {
      stop("'treatment_model' needs to be a model that can be ",
           "refit or a formula object when using bootstrap=TRUE.")
    } else if (is.numeric(obj$outcome_model)) {
      stop("'outcome_model' needs to be an actual model that can be ",
           "refit when using bootstrap=TRUE.")
    }

  }

  # asymptotic variance calculations
  if (conf_int) {
    if (method %in% c("matching", "direct_pseudo")) {
      warning("Asymptotic or exact variance calculations are currently",
              " not available for method='", method, "'. Use bootstrap=TRUE",
              "to get bootstrap estimates.")
    }
  }
}

## for sim_confounded_crisk
check_inputs_sim_crisk_fun <- function(n, lcovars, outcome_betas, gamma,
                                       lambda, treatment_betas, group_beta,
                                       intercept, gtol, cens_fun, cens_args,
                                       max_t, max_iter) {

  if (!is.numeric(n)) {
    stop("'n' must be a positive integer.")
  } else if (floor(n)!=n) {
    stop("'n' must be a positive integer, not a double.")
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
  } else if (!is.function(cens_fun) & !is.null(cens_fun)) {
    stop("'cens_fun' must be a function with the argument 'n' or NULL.")
  } else if (!is.list(cens_args) & !is.null(cens_args)) {
    stop("'cens_args' must be a named list of arguments to be passed to",
         " 'cens_fun' or NULL.")
  } else if (!is.numeric(max_t)) {
    stop("'max_t' must be a number.")
  } else if (max_t <= 0) {
    stop("'max_t' must be bigger than zero.")
  } else if (!is.numeric(group_beta)) {
    stop("'group_beta' must be a number.")
  }

  # TODO: more specific checks here

}
