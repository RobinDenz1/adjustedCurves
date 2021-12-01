## for adjustedsurv function
check_inputs_adjustedsurv <- function(data, variable, ev_time, event, method,
                                      conf_int, conf_level, times, bootstrap,
                                      n_boot, na.action, ...) {
  obj <- list(...)

  if (!inherits(data, c("data.frame", "mids"))) {
    stop("'data' argument must be a data.frame or mids object.")
  # needed variables
  } else if (!is.character(variable) | !is.character(ev_time) |
             !is.character(event) | !is.character(method)) {
    stop("Arguments 'variable', 'ev_time', 'event' and 'method' must be ",
         "character strings, specifying variables in 'data'.")
  # method
  } else if (!method %in% c("km", "iptw_km", "iptw_cox", "iptw_pseudo",
                            "direct", "direct_pseudo", "aiptw_pseudo",
                            "aiptw", "tmle", "ostmle", "matching",
                            "emp_lik", "strat_cupples", "strat_amato",
                            "strat_gregory_nieto")) {
    stop("Method '", method, "' is undefined. See documentation for ",
         "details on available methods.")
  # conf_int
  } else if (!is.logical(conf_int)) {
    stop("'conf_int' must be either TRUE or FALSE.")
  # conf_level
  } else if (!is.numeric(conf_level)) {
    stop("'conf_level' must be a number < 1 and > 0.")
  } else if (conf_level >= 1 | conf_level <= 0) {
    stop("'conf_level' must be a number < 1 and > 0.")
  # bootstrap
  } else if (!is.logical(bootstrap)) {
    stop("'bootstrap' must be either TRUE or FALSE.")
  # n_boot
  } else if (!is.numeric(n_boot) || n_boot < 1) {
    stop("'n_boot' must be a positive integer > 2.")
  # times
  } else if (!is.numeric(times) & !is.null(times)) {
    stop("'times' must be a numeric vector or NULL.")
  }

  # data.frame properties
  if (inherits(data, "data.frame")) {
    # variable
    if (!variable %in% colnames(data)) {
      stop(variable, " is not a valid column name in 'data'.")
    # ev_time
    } else if (!ev_time %in% colnames(data)) {
      stop(ev_time, " is not a valid column name in 'data'.")
    # event
    } else if (!event %in% colnames(data)) {
      stop(event, " is not a valid column name in 'data'.")
    }

    # Check if categorical should be allowed
    levs_len <- length(unique(data[,variable]))
    if (levs_len < 2) {
      stop("There have to be at least two groups in 'variable'.")
    } else if (levs_len > 2 & method %in% c("matching", "emp_lik", "tmle",
                                            "ostmle", "aiptw")) {
      stop("Categorical treatments are currently not supported for ",
           "method='", method, "'.")
    }

    # Check if the group variable has the right format
    if (method %in% c("matching", "emp_lik", "tmle", "ostmle") &
        !is.factor(data[,variable]) & !is.numeric(data[,variable])) {
      stop("The column in 'data' specified by 'variable' needs to be ",
           "a factor or a dichotomous integer variable if method='",
           method, "'.")
    }

    if (!method %in% c("matching", "emp_lik", "tmle", "ostmle") &
        !is.factor(data[,variable])) {
      stop("The column in 'data' specified by 'variable' needs to be ",
           "a factor variable if method='", method, "'.")
    }

    # Check if ev_time variable has the right format
    if (!is.numeric(data[,ev_time])) {
      stop("The column in 'data' specified by 'ev_time' must be numeric.")
    }

    # Check if event variable has the right format
    if (!is.numeric(data[,event]) & !is.logical(data[,event])) {
      stop("The column in 'data' specified by 'event' must be numeric",
           " or logical.")
    }

    # No extrapolation
    if (!is.null(times) && (max(times) > max(data[,ev_time]))) {
      stop("Values in '", ev_time,
           "' must be smaller than max(data[,ev_time]).",
           " No extrapolation allowed.")
    }
  # Check if mira objects where supplied when class(data) == "mids"
  } else {
    # treatment_model
    if (!is.null(obj$treatment_model) &
        !inherits(obj$treatment_model, c("mira", "formula"))) {
      stop("When using multiple imputation, mira objects or a formula",
           " need to be supplied to 'treatment_model' instead of",
           " single models. See documentation.")
    }
    # outcome_model
    if (!is.null(obj$outcome_model) &
        !inherits(obj$outcome_model, "mira")) {
      stop("When using multiple imputation, mira objects need to be supplied",
           " to 'outcome_model' instead of single models. See documentation.")
    }
    # censoring_model
    if (!is.null(obj$censoring_model) &
        !inherits(obj$censoring_model, "mira")) {
      stop("When using multiple imputation, mira objects need to be supplied",
           " to 'censoring_model' instead of single models. See documentation.")
    }
    # warn user when there are missing values in event variable
    if (anyNA(as.data.frame(data$data)[,event])) {
      warning("Using multiple imputation with missing values in 'event'",
              " variable has not been tested yet. Use with caution.",
              call.=FALSE)
    }
    # warn user when there are missing values in ev_time variable
    if (anyNA(as.data.frame(data$data)[,ev_time])) {
      warning("Using multiple imputation with missing values in 'ev_time'",
              " variable has not been tested yet. Use with caution.",
              call.=FALSE)
    }
  }

  ## Direct Pseudo, AIPTW Pseudo
  if (method=="direct_pseudo" | method=="aiptw_pseudo" |
      method=="iptw_pseudo") {

    # outcome_vars
    if ("outcome_vars" %in% names(obj) && (!is.character(obj$outcome_vars))) {
      stop("'outcome_vars' should be a character vector of column names ",
           "in 'data', used to model the outcome mechanism.")
    }
    # type_time
    if ("type_time" %in% names(obj) && (!obj$type_time %in%
                                        c("factor", "bs", "ns"))) {
      stop("'type_time' should be either 'factor', 'bs' or 'ns'.")
    }
    # treatment_model
    if ("treatment_model" %in% names(obj) && method=="aiptw_pseudo" &
        (!(is.numeric(obj$treatment_model) |
           inherits(obj$treatment_model, c("glm", "multinom", "mira"))))) {
      stop("Argument 'treatment_model' must be one of: 'glm',",
           " 'multinom', 'mira' or a numeric vector of propensity scores.")
    }

    if (method %in% c("aiptw_pseudo", "direct_pseudo") & !is.null(times)) {
      # times
      if (length(times)==1) {
        stop("'geese' models require at least two distinct time points. ",
             "Add more points in time to 'times' and try again.")
      # spline_df
      } else if (!is.null(obj$spline_df) && (obj$spline_df > length(times))) {
        warning("'spline_df' > len(times) might lead to problems.",
                call.=FALSE)
      } else if (!is.null(obj$type_time) && (5 > length(times) &
                                             obj$type_time!="factor")) {
        warning("'spline_df' > length(times)=5 might lead to problems when",
                  " type_time!='factor'.", call.=FALSE)
      }
    }

  ## TMLE / OSTMLE
  } else if (method=="tmle" | method=="ostmle") {
    # times
    if (!is.null(times) && (!all(times==floor(times)))) {
      stop("Only integer time is allowed when using method='tmle' or",
           " method='ostmle'.")
    }
  ## Empirical Likelihood
  } else if (method=="emp_lik") {

    # treatment_vars
    if (!is.character(obj$treatment_vars)) {
      stop("'treatment_vars' should be a character vector of column names ",
           "in 'data', used to model the outcome mechanism.")
    # moment
    } else if ("moment" %in% names(obj) && (!obj$moment %in%
                                            c("first", "second"))) {
      stop("Argument 'moment' must be either 'first' or 'second'.")
    # standardize
    } else if ("standardize" %in% names(obj) &&
               (!is.logical(obj$standardize))) {
      stop("Argument 'standardize' must be either TRUE or FALSE.")
    } else if (inherits(data, "data.frame")) {

      # treatment_vars
      if (!all(obj$treatment_vars %in% colnames(data))) {
        stop("'treatment_vars' should be a character vector of column names ",
             "in 'data', used to model the outcome mechanism.")
      }

      for (col in obj$treatment_vars) {
        if (paste0(unique(data[,col]), collapse="") %in% c("10", "01")) {
          warning("Dichotomous variables coded with 0 and 1 found in ",
                  " 'treatment_vars'. Consider recoding to -1 and 1",
                  " to avoid estimation problems.", call.=FALSE)
        }
      }

    }
  ## Matching
  } else if (method=="matching") {

    # treatment model
    if (!"treatment_model" %in% names(obj)) {
      stop("Argument 'treatment_model' must be specified when using",
           " method='matching'.")
    } else if (is.numeric(obj$treatment_model) &&
               (any(obj$treatment_model > 1) | any(obj$treatment_model < 0))) {
      stop("Propensity Scores > 1 or < 0 not allowed. Perhaps you supplied ",
           "weights on accident?")
    }
  ## AIPTW
  } else if (method=="aiptw") {
    # need treatment_model
    if (!"treatment_model" %in% names(obj)) {
      stop("Argument 'treatment_model' must be specified when using",
           " method='aiptw'.")
    # need outcome_model
    } else if (!"outcome_model" %in% names(obj)) {
      stop("Argument 'outcome_model' must be specified when using",
           " method='aiptw'.")
    }
  ## Direct
  } else if (method=="direct") {
    # need outcome_model
    if (!"outcome_model" %in% names(obj)) {
      stop("Argument 'outcome_model' must be specified when using",
           " method='direct'.")
    # no bootstrapping with pecRpart
    } else if (inherits(obj$outcome_model, c("pecRpart")) & bootstrap) {
      stop("Bootstrapping is currently not supported with method='direct'",
           " and an 'outcome_model' of class 'pecRpart'.")
    # only allow certain models when there is no censoring
    } else if (inherits(obj$outcome_model, c("glm", "ols", "randomForest")) &&
               all(data[,event]==1)) {
      stop("'outcome_models' of class c('glm', 'ols', 'randomForest') are",
           " only allowed when there is no censoring.")
    # don't allow selectCox if no covariates are left after selection
    } else if (inherits(obj$outcome_model, c("selectCox")) &&
               inherits(obj$outcome_model$fit, "survfit")) {
      stop("The final 'fit' object in the 'selectCox' model must be a",
           " coxph object, not survfit.")
    # variable has to be inside formula of coxph
    } else if (inherits(obj$outcome_model, "coxph")) {
      if (!variable %in% all.vars(obj$outcome_model$formula)) {
        stop("'variable' has to be included in the cox-regression model.")
      }
    }
  ## Cupples / Amato / Gregory
  } else if (method=="strat_cupples" | method=="strat_amato" |
             method=="strat_gregory_nieto") {
    # need adjust_vars
    if (!"adjust_vars" %in% names(obj)) {
      stop("Argument 'adjust_vars' needs to be specified when using",
           " method='", method, "'.")
    # no continuous confounders
    } else if (inherits(data, "data.frame")) {
      for (i in seq_len(length(obj$adjust_vars))) {
        if (is.numeric(data[,obj$adjust_vars[i]]) &&
            !all(floor(data[,obj$adjust_vars[i]])==data[,obj$adjust_vars[i]])) {
          stop("Variables in 'adjust_vars' have to be integer, factor or",
               "character variables. Continuous variables are not allowed",
               " when using method='", method, "'.")
        }
      }
    }
    # valid reference data
    if ((method=="strat_cupples" | method=="strat_amato") &
        !is.null(obj$reference)) {
      if (!all(obj$adjust_vars %in% colnames(obj$reference))) {
        stop("If a 'reference' data.frame is supplied, it needs to contain",
             " all variables listed in 'adjust_vars'.")
      }
    }
  }

  # bootstrapping
  if (bootstrap & is.numeric(obj$treatment_model)) {
    stop("'treatment_model' needs to be a model that can be ",
         "refit or a formula object when using bootstrap=TRUE.")
  }

  # asymptotic variance calculations
  if (conf_int & (method %in% c("emp_lik", "matching", "direct_pseudo",
                                "strat_cupples", "strat_amato"))) {
    warning("Asymptotic or exact variance calculations are currently",
            " not available for method='", method, "'. Use bootstrap=TRUE",
            " to get bootstrap estimates.", call.=FALSE)
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
        (is.list(lcovars) & is.numeric(outcome_betas) &
         is.numeric(treatment_betas)))) {
    stop("'lcovars', 'outcome_betas' and 'treatment_betas' can either all ",
         "be NULL to use default values, or have to be specified.")
  }

  if (!is.null(lcovars) && (is.null(names(lcovars)))) {
    stop("Elements in the 'lcovars' list must be named.")
  }

  if (!is.null(outcome_betas) && (is.null(names(outcome_betas)))) {
    stop("Elements in the 'outcome_betas' vector must be named.")
  }

  if (!is.null(treatment_betas) && (is.null(names(treatment_betas)))) {
    stop("Elements in the 'treatment_betas' vector must be named.")
  }

  lens <- c(length(outcome_betas), length(treatment_betas), length(lcovars))
  if (!is.null(lcovars) && stats::var(lens)!=0) {
    stop("'outcome_betas', 'treatment_betas' and 'lcovars' must have",
         " the same length.")
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
    for (i in seq_len(length(lcovars))) {
      if (!lcovars[[i]][1] %in% c("rbinom", "rnorm", "runif")) {
        stop("The first element of every vector in 'lcovars' must be either ",
             "'rbinom', 'rnorm' or 'runif', not ", lcovars[[i]][1], ".")
      }
    }
  }
}

## for adjusted_rmst function
check_inputs_adj_rmst <- function(adjsurv, from, to, use_boot) {

  if (!is.numeric(from) | !is.numeric(to)) {
    stop("'from' and 'to' must be numbers.")
  } else if (!(from >= 0 & to >= 0)) {
    stop("'from' and 'to' must be >= 0.")
  } else if (!inherits(adjsurv, "adjustedsurv")) {
    stop("'adjsurv' must be an 'adjustedsurv' object created using ",
         "the 'adjustedsurv()' function.")
  } else if (from >= to) {
    stop("'from' must be smaller than 'to'.")
  } else if (use_boot & is.null(adjsurv$boot_adjsurv)) {
    warning("Cannot use bootstrapped estimates because",
            " they were not estimated.",
            " Need 'bootstrap=TRUE' in 'adjustedsurv' function call.",
            call.=FALSE)
  }

  if (to > max(adjsurv$adjsurv$time, na.rm=TRUE)) {
    stop("'to' can not be greater than the latest observed time.")
  }

  if (length(unique((adjsurv$adjsurv$time))) < 10) {
    warning("Using only a few points in time might lead to biased",
            " estimates. Consider using a finer times grid in",
            " 'adjustedsurv'.", call.=FALSE)
  }

}

## for adjusted_rmtl function
check_inputs_adj_rmtl <- function(adj, from, to, use_boot) {

  if (!is.numeric(from) | !is.numeric(to)) {
    stop("'from' and 'to' must be numbers.")
  } else if (!(from >= 0 & to >= 0)) {
    stop("'from' and 'to' must be >= 0.")
  } else if (!inherits(adj, c("adjustedsurv", "adjustedcif"))) {
    stop("'adj' must be an 'adjustedsurv' object created using",
         " the 'adjustedsurv()' function or an 'adjustedcif' object",
         " created using the 'adjustedcif()' function.")
  } else if (from >= to) {
    stop("'from' must be smaller than 'to'.")
  } else if (use_boot & is.null(adj$boot_adjsurv) & is.null(adj$boot_adjcif)) {
    warning("Cannot use bootstrapped estimates because",
            " they were not estimated.",
            " Need 'bootstrap=TRUE' in 'adjustedsurv'/'adjustedcif'",
            " function call.", call.=FALSE)
  }

  if (inherits(adj, "adjustedsurv")) {
    max_t <- max(adj$adjsurv$time, na.rm=TRUE)
    n_t <- length(unique((adj$adjsurv$time)))
  } else {
    max_t <- max(adj$adjcif$time, na.rm=TRUE)
    n_t <- length(unique((adj$adjcif$time)))
  }

  if (to > max_t) {
    stop("'to' can not be greater than the latest observed time.")
  }

  if (n_t < 10) {
    warning("Using only a few points in time might lead to biased",
            " estimates. Consider using a finer times grid in",
            " 'adjustedsurv'/'adjustedcif'.", call.=FALSE)
  }
}

## for adjustedsurv_test
check_inputs_adj_test <- function(adjsurv, from, to) {

  if (!(inherits(adjsurv, "adjustedsurv") |
        inherits(adjsurv, "adjustedcif"))) {
    stop("'adjsurv' must be an 'adjustedsurv' or 'adjustedcif' object,",
         "created using the adjustedsurv or adjustedcif function.")
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

  if (inherits(adjsurv, "adjustedsurv")) {
    if (to > max(adjsurv$adjsurv$time, na.rm=TRUE)) {
      stop("'to' can not be greater than the latest observed time.")
    }
  } else {
    if (to > max(adjsurv$adjcif$time)) {
      stop("'to' can not be greater than the latest observed time.")
    }
  }

  if (inherits(adjsurv, "adjustedsurv")) {
    if (length(unique((adjsurv$adjsurv$time))) < 10) {
      warning("Using only a few points in time might lead to biased",
              " estimates. Consider using a finer times grid in",
              " 'adjustedsurv'.", call.=FALSE)
    }
  } else {
    if (length(unique((adjsurv$adjcif$time))) < 10) {
      warning("Using only a few points in time might lead to biased",
              " estimates. Consider using a finer times grid in",
              " 'adjustedcif'.", call.=FALSE)
    }
  }
}

## for adjustedcif
check_inputs_adjustedcif <- function(data, variable, ev_time, event, method,
                                     conf_int, conf_level, times, bootstrap,
                                     n_boot, cause=cause, na.action, ...) {
  obj <- list(...)

  if (!inherits(data, c("data.frame", "mids"))) {
    stop("'data' argument must be a data.frame or mids object.")
  # needed variables
  } else if (!is.character(variable) | !is.character(ev_time) |
             !is.character(event) | !is.character(method)) {
    stop("Arguments 'variable', 'ev_time', 'event' and 'method' must be ",
         "character strings, specifying variables in 'data'.")
  } else if (!method %in% c("aalen_johansen", "iptw", "iptw_pseudo", "direct",
                            "direct_pseudo", "aiptw_pseudo",
                            "aiptw", "tmle", "matching")) {
    stop("Method '", method, "' is undefined. See documentation for ",
         "details on available methods.")
  # conf_int
  } else if (!is.logical(conf_int)) {
    stop("'conf_int' must be either TRUE or FALSE.")
  # conf_level
  } else if (!is.numeric(conf_level)) {
    stop("'conf_level' must be a number < 1 and > 0.")
  } else if (conf_level >= 1 | conf_level <= 0) {
    stop("'conf_level' must be a number < 1 and > 0.")
  # cause
  } else if (!is.numeric(cause)) {
    stop("'cause' must be a number specifying the cause of interest in ",
         "the column specified with 'event'.")
  } else if (length(cause)!=1) {
    stop("'cause' must be of length = 1.")
  # time
  } else if (!is.numeric(times) & !is.null(times)) {
    stop("'times' must be a numeric vector or NULL.")
  # bootstrap
  } else if (!is.logical(bootstrap)) {
    stop("'bootstrap' must be either TRUE or FALSE.")
  # n_boot
  } else if (!is.numeric(n_boot) || n_boot < 2) {
    stop("'n_boot' must be a positive integer > 1.")
  }

  # data.frame properties
  if (inherits(data, "data.frame")) {
    # variable
    if (!variable %in% colnames(data)) {
      stop(variable, " is not a valid column name in 'data'.")
    # ev_time
    } else if (!ev_time %in% colnames(data)) {
      stop(ev_time, " is not a valid column name in 'data'.")
    # event
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

    # Check if ev_time variable has the right format
    if (!is.numeric(data[,ev_time])) {
      stop("The column in 'data' specified by 'ev_time' must be numeric.")
    }

    # Check if event variable has the right format
    if (!is.numeric(data[,event]) & !is.logical(data[,event])) {
      stop("The column in 'data' specified by 'event' must be numeric",
           " or logical.")
    }

    # Check if categorical should be allowed
    levs_len <- length(unique(data[,variable]))
    if (levs_len < 2) {
      stop("There have to be at least two groups in 'variable'.")
    } else if (levs_len > 2 & method %in% c("matching", "tmle", "aiptw")) {
      stop("Categorical treatments are currently not supported for ",
           "method='", method, "'.")
    }

    # No extrapolation
    if (!is.null(times) && (max(times) > max(data[,ev_time]))) {
      stop("Values in '", ev_time,
           "' must be smaller than max(data[,ev_time]).",
           " No extrapolation allowed.")
    }
  }

  ## Direct Pseudo, AIPTW Pseudo
  if (method=="direct_pseudo" | method=="aiptw_pseudo" |
      method=="iptw_pseudo") {

    # outcome_vars
    if ("outcome_vars" %in% names(obj) && (!is.character(obj$outcome_vars))) {
      stop("'outcome_vars' should be a character vector of column names",
           " in 'data', used to model the outcome mechanism.")
    }
    # type_time
    if ("type_time" %in% names(obj) && (!obj$type_time %in%
                                        c("factor", "bs", "ns"))) {
      stop("'type_time' should be either 'factor', 'bs' or 'ns'.")
    }
    # treatment_model
    if ("treatment_model" %in% names(obj) && method=="aiptw_pseudo" &
        (!(is.numeric(obj$treatment_model) |
           inherits(obj$treatment_model, c("glm", "multinom", "mira"))))) {
      stop("Argument 'treatment_model' must be one of",
           " glm, multinom or mira or a numeric vector of propensity scores.")
    }

    if (method %in% c("aiptw_pseudo", "direct_pseudo") & !is.null(times)) {
      # times
      if (length(times)==1) {
        stop("'geese' models require at least two distinct time points. ",
             "Add more points in time to 'times' and run again.")
      # spline_df
      } else if (!is.null(obj$spline_df) && (obj$spline_df > length(times))) {
        warning("'spline_df' > len(times) might lead to problems.",
                call.=FALSE)
      } else if (!is.null(obj$type_time) &&
                 (5 > length(times) & obj$type_time!="factor")) {
        warning("'spline_df' > len(times) might lead to problems.",
                call.=FALSE)
      }
    }

  ## TMLE
  } else if (method=="tmle") {
    if (!is.null(times) && (!all(times==floor(times)))) {
      stop("Only integer time is allowed when using method='tmle'.")
    }
  ## Matching
  } else if (method=="matching") {
    # treatment_model
    if (!"treatment_model" %in% names(obj)) {
      stop("Argument 'treatment_model' must be specified when using",
           " method='matching'.")
    } else if (is.numeric(obj$treatment_model) &&
               (any(obj$treatment_model > 1) | any(obj$treatment_model < 0))) {
      stop("Propensity Scores > 1 or < 0 not allowed. Perhaps you supplied ",
           "weights on accident?")
    }
  ## AIPTW
  } else if (method=="aiptw") {
    # need treatment_model
    if (!"treatment_model" %in% names(obj)) {
      stop("Argument 'treatment_model' must be specified when using",
           " method='aiptw'.")
      # need outcome_model
    } else if (!"outcome_model" %in% names(obj)) {
      stop("Argument 'outcome_model' must be specified when using",
           " method='aiptw'.")
    }
  ## Direct
  } else if (method=="direct") {
    if (!"outcome_model" %in% names(obj)) {
      stop("Argument 'outcome_model' must be specified when using",
           " method='direct'.")
    }
    if (inherits(obj$outcome_model, "FGR") &&
        cause != obj$outcome_model$cause) {
      stop("The FGR model needs to be fit with the same 'cause' as specified",
           " in the 'cause' argument.")
    }
  }

  # bootstrapping
  if (bootstrap && is.numeric(obj$treatment_model)) {
    stop("'treatment_model' needs to be a model that can be",
         " refit or a formula object when using bootstrap=TRUE.")
  }

  # asymptotic variance calculations
  if (conf_int & (method %in% c("matching", "direct_pseudo"))) {
    warning("Asymptotic or exact variance calculations are currently",
            " not available for method='", method, "'. Use bootstrap=TRUE",
            " to get bootstrap estimates.", call.=FALSE)
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
    stop("'group_beta' must be a numeric vector.")
  }

  if (stats::var(c(length(group_beta), length(gamma), length(lambda)))!=0) {
    stop("Arguments 'group_beta', 'gamma' and 'lambda' need to have",
         " the same length.")
  }

  if (!((is.null(lcovars) & is.null(outcome_betas) & is.null(treatment_betas)) |
        (is.list(lcovars) & is.list(outcome_betas) &
         is.numeric(treatment_betas)))) {
    stop("'lcovars', 'outcome_betas' and 'treatment_betas' can either all ",
         "be NULL to use default values, or have to be specified.")
  }

  lens <- c(length(outcome_betas), length(treatment_betas), length(lcovars))
  if (!is.null(lcovars) && stats::var(lens)!=0) {
    stop("'outcome_betas', 'treatment_betas' and 'lcovars' must have",
         " the same length.")
  }

  if (!is.null(outcome_betas) && stats::var(vapply(outcome_betas, length,
                                                   FUN.VALUE=numeric(1)))!=0) {
    stop("The vectors supplied in 'outcome_betas' all need to have,",
         " the same length.")
  }

  if (!is.null(lcovars) && (is.null(names(lcovars)))) {
    stop("Elements in the 'lcovars' list must be named.")
  }

  if (!is.null(treatment_betas) && (is.null(names(treatment_betas)))) {
    stop("Elements in the 'treatment_betas' vector must be named.")
  }

  if (!is.null(lcovars)) {
    for (i in seq_len(length(lcovars))) {
      if (!lcovars[[i]][1] %in% c("rbinom", "rnorm", "runif")) {
        stop("The first element of every vector in 'lcovars' must be either ",
             "'rbinom', 'rnorm' or 'runif', not ", lcovars[[i]][1], ".")
      }
    }
  }
}
