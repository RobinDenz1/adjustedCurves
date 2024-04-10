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

## estimate iptw weights
get_iptw_weights <- function(data, treatment_model, weight_method,
                             variable, stabilize=TRUE, trim, trim_q, ...) {
  levs <- levels(data[, variable])

  # using WeightIt
  if (inherits(treatment_model, "formula")) {
    args <- list(formula=treatment_model, data=data,
                 method=weight_method,
                 estimand="ATE")
    weights <- do.call(WeightIt::weightit, c(args, ...))$weights
  # using a logistic regression model
  } else if (inherits(treatment_model, "glm")) {
    ps <- stats::predict.glm(treatment_model, newdata=data, type="response")
    weights <- ifelse(data[, variable]==levs[2], 1/ps, 1/(1-ps))
  # using a multinomial logistic regression model
  } else if (inherits(treatment_model, "multinom")) {
    predict.multinom <- utils::getFromNamespace("predict.multinom", "nnet")
    ps <- predict.multinom(treatment_model, newdata=data, type="probs")

    weights <- rep(0, nrow(data))
    for (i in levs) {
      weights[data[, variable] == i] <- 1/ps[data[, variable] == i, i]
    }
  # user supplied a mira model but no mids data
  } else if (inherits(treatment_model, "mira")) {
    stop("If a 'mira' object is used in the 'treatment_model' argument",
         " the 'data' argument must be a 'mids' object, not a",
         " data.frame.", call.=FALSE)
  # nothing else allowed
  } else {
    stop("Unsuported input: '", class(treatment_model),
         "'. See documentation.")
  }

  weights <- trim_weights(weights=weights, trim=trim)
  weights <- trim_weights_quantiles(weights=weights, trim_q=trim_q)

  if (stabilize) {
    weights <- stabilize_weights(weights, data, variable, levs)
  }

  return(weights)
}

## trim weights that are above a certain value
trim_weights <- function(weights, trim) {

  if (!trim) {
    return(weights)
  } else {
    weights[weights > trim] <- trim
    return(weights)
  }
}

## trim weights based on defined quantiles
trim_weights_quantiles <- function(weights, trim_q) {

  # check inputs
  if (!((length(trim_q)==1 && is.logical(trim_q) &&
         !trim_q) | (length(trim_q)==2 &&
         is.numeric(trim_q) && all(trim_q < 1)
         && all(trim_q > 0)))) {
    stop("'trim_quantile' must be either FALSE or a numeric vector of length 2",
         " containing the lower and upper quantile to be trimmed.")
  }

  if (length(trim_q)==1) {
    return(weights)
  } else {
    q_low <- stats::quantile(weights, probs=min(trim_q))
    q_high <- stats::quantile(weights, probs=max(trim_q))
    weights[weights < q_low] <- q_low
    weights[weights > q_high] <- q_high
    return(weights)
  }
}

## stabilize weights
stabilize_weights <- function(weights, data, variable, levs) {

  w_names <- names(weights)
  tab <- stats::setNames(vapply(X=levs, FUN=function(x, treat) {mean(treat==x)},
                                FUN.VALUE=numeric(1L), treat=data[,variable]),
                         levs)
  weights <- stats::setNames(weights * tab[as.character(data[,variable])],
                             w_names)

  return(weights)
}

## function to get predicted values from any kind of geese object,
## for multiple time points
geese_predictions <- function(geese_mod, Sdata, times, n) {

  # full model matrix and betas
  mod_mat <- stats::model.matrix(geese_mod$formula,
                                 data=stats::model.frame(geese_mod$formula,
                                             Sdata,
                                             na.action=stats::na.pass))
  betas <- geese_mod$beta

  apply_betas <- function(x, betas) {
    return(sum(x * betas))
  }

  pred_mat <- matrix(nrow=n, ncol=length(times))
  for (i in seq_len(length(times))) {
    # take only relevant portion (at time t) of model matrix
    mod_mat_t <- mod_mat[Sdata$vtime==times[i], ]
    # apply coefficients
    preds <- apply(X=mod_mat_t, MARGIN=1, FUN=apply_betas, betas=betas)
    pred_mat[, i] <- preds
  }
  colnames(pred_mat) <- paste0("t.", times)
  return(pred_mat)
}

## calculate CI from sd
confint_surv <- function(surv, se, conf_level, conf_type="plain") {
  # critical value
  z_val <- stats::qnorm(1-((1-conf_level)/2))

  if (conf_type=="plain") {
    error <- z_val * se
    left <- surv - error
    right <- surv + error
  } else if (conf_type=="log") {
    error <- z_val * se
    left <- surv * exp(-error)
    right <- surv * exp(error)
  }

  return(list(left=left, right=right))
}

## function to change plotdata in iptw and standard methods when
## custom points in time are supplied
specific_times <- function(plotdata, times, est="surv", interpolation="steps") {

  levs <- unique(plotdata$group)
  new_plotdata <- vector(mode="list", length=length(levs))
  for (i in seq_len(length(levs))) {

    new_est <- read_from_fun(x=times, est=est,
                             data=plotdata[which(plotdata$group==levs[i]), ],
                             interpolation=interpolation)

    if (est=="cif") {
      new_dat <- data.frame(time=times, cif=new_est, group=levs[i])
    } else if (est=="surv") {
      new_dat <- data.frame(time=times, surv=new_est, group=levs[i])
    } else if (est=="diff") {
      new_dat <- data.frame(time=times, diff=new_est, group=levs[i])
    } else if (est=="ratio") {
      new_dat <- data.frame(time=times, ratio=new_est, group=levs[i])
    }

    if ("se" %in% colnames(plotdata)) {
      # read from curve using custom function
      # NOTE: while the confidence intervals should be interpolated according
      #       to the respective specification, the standard error stays
      #       the same over the interpolation period and therefore always
      #       needs to be interpolated using the step function method
      new_se <- read_from_fun(x=times, interpolation="steps", est="se",
                       data=plotdata[which(plotdata$group==levs[i]), ])
      new_ci_lower <- read_from_fun(x=times, interpolation=interpolation,
                                    est="ci_lower",
                              data=plotdata[which(plotdata$group==levs[i]), ])
      new_ci_upper <- read_from_fun(x=times, interpolation=interpolation,
                                    est="ci_upper",
                              data=plotdata[which(plotdata$group==levs[i]), ])
      # add to output in same order
      new_dat$se <- new_se
      new_dat$ci_lower <- new_ci_lower
      new_dat$ci_upper <- new_ci_upper
    }

    # same as standard error, stays the same
    if ("p_value" %in% colnames(plotdata)) {
      new_p <- read_from_fun(x=times, interpolation="steps", est="p_value",
                      data=plotdata[which(plotdata$group==levs[i]), ])
      new_dat$p_value <- new_p
    }

    new_plotdata[[i]] <- new_dat
  }
  new_plotdata <- as.data.frame(dplyr::bind_rows(new_plotdata))
  return(new_plotdata)
}

## calculate pseudo values for the survival function
## using either the standard way or with dependent censoring
calc_pseudo_surv <- function(data, ev_time, event, times, censoring_vars,
                             ipcw.method, cause=1) {

  # standard pseudo-values, no dependent censoring
  if (is.null(censoring_vars)) {

    # estimate pseudo observations
    hist_formula <- stats::as.formula(paste("prodlim::Hist(", ev_time, ", ",
                                            event, ") ~ 1"))
    pseudo <- prodlim::jackknife(prodlim::prodlim(hist_formula, data=data),
                                 times=times, cause=cause)

  } else {

    requireNamespace("eventglm")
    pseudo_aareg <- utils::getFromNamespace("pseudo_aareg", "eventglm")

    cens_formula <- stats::as.formula(paste0("~ ", paste(censoring_vars,
                                                         collapse=" + ")))
    pseudo_formula <- stats::as.formula(paste0("survival::Surv(", ev_time,
                                               ", ", event, ") ~ 1"))

    pseudo <- vapply(times, FUN=pseudo_aareg,
                     formula=pseudo_formula,
                     cause=1,
                     data=data,
                     type="survival",
                     formula.censoring=cens_formula,
                     ipcw.method=ipcw.method,
                     FUN.VALUE=numeric(nrow(data)))
  }
  return(pseudo)
}

## keep only the covariates needed for the analysis
## this has to be done in order to correctly use na.action
remove_unnecessary_covars <- function(data, method, variable, ev_time,
                                      event, ...) {

  # nothing is removed for tmle
  if (method=="tmle") {
    return(data)
  }

  args <- list(...)

  # extract variables from treatment model
  if (inherits(args$treatment_model, "multinom")) {
    treatment_vars <- all.vars(args$treatment_model$call$formula)
  } else if (inherits(args$treatment_model, c("glm", "lm"))) {
    treatment_vars <- all.vars(args$treatment_model$formula)
  } else if (inherits(args$treatment_model, "formula")) {
    treatment_vars <- all.vars(args$treatment_model)
  } else {
    treatment_vars <- NULL
  }

  # extract variables from outcome model
  if (inherits(args$outcome_model, c("coxph", "mexhaz"))) {
    outcome_vars <- all.vars(args$outcome_model$formula)
  } else if (inherits(args$outcome_model, c("CauseSpecificCox", "FGR", "aalen",
                                            "cox.aalen", "flexsurvreg",
                                            "pecCforest", "prodlim",
                                            "psm", "randomForest",
                                            "riskRegression", "selectCox",
                                            "glm", "ols", "rfsrc",
                                            "penfitS3", "gbm",
                                            "singleEventCB", "fcrr",
                                            "comprisk"))) {
    outcome_vars <- all.vars(args$outcome_model$call$formula)
  } else if (inherits(args$outcome_model, "pecRpart")) {
    outcome_vars <- all.vars(args$outcome_model$rpart$terms)
  } else if (inherits(args$outcome_model, "ranger")) {
    outcome_vars <- all.vars(args$outcome_model$call[[2]])
  } else {
    outcome_vars <- NULL
  }

  # extract variables from censoring model
  if (inherits(args$censoring_model, "coxph")) {
    censoring_vars <- all.vars(args$censoring_model$formula)
  } else {
    censoring_vars <- NULL
  }

  # covariates that are always needed
  needed_covars <- c(variable, ev_time, event)

  # method specific covariate needs
  if (method=="direct") {
    needed_covars <- c(needed_covars, outcome_vars)
  } else if (method=="direct_pseudo") {
    needed_covars <- c(needed_covars, args$outcome_vars, args$censoring_vars)
  } else if (method %in% c("iptw", "iptw_km", "iptw_cox", "iptw_pseudo",
                           "matching")) {
    needed_covars <- c(needed_covars, treatment_vars, args$censoring_vars)
  } else if (method=="emp_lik") {
    needed_covars <- c(needed_covars, args$treatment_vars)
  } else if (method=="aiptw") {
    needed_covars <- c(needed_covars, treatment_vars, outcome_vars,
                       censoring_vars)
  } else if (method=="aiptw_pseudo") {
    needed_covars <- c(needed_covars, args$outcome_vars, args$censoring_vars,
                       treatment_vars)
  } else if (method=="strat_cupples" | method=="strat_amato" |
             method=="strat_nieto") {
    needed_covars <- c(needed_covars, args$adjust_vars)
  } else if (method=="iv_2SRIF") {
    needed_covars <- c(needed_covars, args$adjust_vars, args$instrument)
  } else if (method=="prox_iptw" | method=="prox_aiptw") {
    needed_covars <- c(needed_covars, args$adjust_vars, args$treatment_proxy,
                       args$outcome_proxy)
  }

  # remove duplicates
  needed_covars <- unique(needed_covars)

  # filter data
  data <- dplyr::select(data, dplyr::all_of(needed_covars))

  return(data)
}

## require needed packages
load_needed_packages <- function(method, kind, treatment_model,
                                 censoring_vars) {

  if (kind=="surv") {

    # survival
    if (method=="direct" | method=="km" | method=="strat_cupples" |
        method=="tmle" | method=="iv_2SRIF") {
      requireNamespace("survival")
    }

    # riskRegression
    if (method=="aiptw" | method=="direct") {
      requireNamespace("riskRegression")
    }

    # pseudo-values
    if (method %in% c("direct_pseudo", "aiptw_pseudo", "iptw_pseudo")) {
      requireNamespace("prodlim")

      if (!is.null(censoring_vars)) {
        requireNamespace("eventglm")
      }
    }

    # WeightIt
    if ((method %in% c("iptw_km", "iptw_cox", "iptw_pseudo"))
               && inherits(treatment_model, "formula")) {
      requireNamespace("WeightIt")
    }

    # multinom
    if ((method %in% c("iptw_km", "iptw_cox", "iptw_pseudo", "aiptw_pseudo"))
               && inherits(treatment_model, "multinom")) {
      requireNamespace("nnet")
    }

    # geese
    if (method=="direct_pseudo" | method=="aiptw_pseudo") {
      requireNamespace("geepack")
    }

    # MASS
    if (method=="emp_lik") {
      requireNamespace("MASS")
    }

    # concrete
    if (method=="tmle") {
      #requireNamespace("concrete")
    }

    # data.table
    if (method=="tmle") {
      requireNamespace("data.table")
    }

    if (method=="prox_iptw" | method=="prox_aiptw") {
      requireNamespace("numDeriv")
    }

  } else {

    # cmprsk
    if (method=="aalen_johansen") {
      requireNamespace("cmprsk")
    }

    # riskRegression
    if (method=="aiptw" | method=="direct" | method=="iptw") {
      requireNamespace("riskRegression")
    }

    # pseudo-values
    if (method %in% c("direct_pseudo", "aiptw_pseudo", "iptw_pseudo",
                      "direct")) {
      requireNamespace("prodlim")
    }

    # WeightIt
    if (method=="iptw_pseudo" && inherits(treatment_model, "formula")) {
      requireNamespace("WeightIt")
    }

    # multinom
    if ((method %in% c("iptw", "iptw_pseudo", "aiptw_pseudo"))
        && inherits(treatment_model, "multinom")) {
      requireNamespace("nnet")
    }

    # geese
    if (method=="direct_pseudo" | method=="aiptw_pseudo") {
      requireNamespace("geepack")
    }

    # concrete
    if (method=="tmle") {
      #requireNamespace("concrete")
    }

    # data.table
    if (method=="tmle") {
      requireNamespace("data.table")
    }
  }
}

## suppress cat() output
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

## add rows with zero survival time if needed for plot
add_rows_with_zero <- function(plotdata, mode="surv") {

  . <- group <- time <- no_zero <- NULL

  # which groups have no zero?
  no_zero_dat <- plotdata %>%
    dplyr::group_by(., group) %>%
    dplyr::summarise(no_zero=!any(time==0, na.rm=TRUE)) %>%
    dplyr::filter(., no_zero)
  levs_no_zero <- unique(no_zero_dat$group)

  if (length(levs_no_zero)!=0 & mode=="surv") {
    row_0 <- data.frame(time=0, group=levs_no_zero, surv=1)

    if ("ci_lower" %in% colnames(plotdata)) {
      row_0$se <- 0
      row_0$ci_lower <- 1
      row_0$ci_upper <- 1

      if ("boot_surv" %in% colnames(plotdata)) {
        row_0$boot_surv <- 1
        row_0$n_boot <- plotdata$n_boot[1]
        row_0 <- dplyr::select(row_0, c("time", "group", "boot_surv",
                                        "se", "ci_lower", "ci_upper",
                                        "n_boot", "surv"))
      }
    }
    rownames(row_0) <- NULL
    plotdata <- rbind(row_0, plotdata)

  } else if (length(levs_no_zero)!=0 & mode=="cif") {
    row_0 <- data.frame(time=0, group=levs_no_zero, cif=0)

    if ("ci_lower" %in% colnames(plotdata)) {
      row_0$se <- 0
      row_0$ci_lower <- 0
      row_0$ci_upper <- 0

      if ("boot_cif" %in% colnames(plotdata)) {
        row_0$boot_cif <- 0
        row_0$n_boot <- plotdata$n_boot[1]
        row_0 <- dplyr::select(row_0, c("time", "group", "boot_cif",
                                        "se", "ci_lower", "ci_upper",
                                        "n_boot", "cif"))
      }
    }
    rownames(row_0) <- NULL
    plotdata <- rbind(row_0, plotdata)
  }

  return(plotdata)
}

## perform isotonic regression on survival / CIF estimates
iso_reg_est <- function(plotdata) {

  mode <- ifelse("surv" %in% colnames(plotdata), "surv", "cif")

  if (anyNA(plotdata[, mode])) {
    stop("Isotonic Regression cannot be used when there are missing",
         " values in the final estimates.")
  }

  for (lev in levels(plotdata$group)) {
    temp <- plotdata[plotdata$group==lev, ]

    if (mode=="surv") {
      new <- rev(stats::isoreg(rev(temp$surv))$yf)
    } else {
      new <- stats::isoreg(temp$cif)$yf
    }
    plotdata[, mode][plotdata$group==lev] <- new

    # shift confidence intervals accordingly
    if ("ci_lower" %in% colnames(temp)) {
      diff <- temp[, mode] - new

      plotdata$ci_lower[plotdata$group==lev] <- temp$ci_lower - diff
      plotdata$ci_upper[plotdata$group==lev] <- temp$ci_upper - diff
    }
  }

  return(plotdata)
}

## force probabilities to be in the 0/1 range
force_bounds_est <- function(plotdata) {
  mode <- ifelse("surv" %in% colnames(plotdata), "surv", "cif")

  if (mode=="surv") {
    plotdata <- within(plotdata, {
      surv <- ifelse(surv < 0, 0, surv)
      surv <- ifelse(surv > 1, 1, surv)
    })
  } else {
    plotdata <- within(plotdata, {
      cif <- ifelse(cif < 0, 0, cif)
      cif <- ifelse(cif > 1, 1, cif)
    })
  }

  return(plotdata)
}

## Given a mids object and our column names of interest, calculate the
## maximum observed (cause-specific) event time
max_observed_time <- function(mids, variable, ev_time, event, levs, cause,
                              method, type) {

  if (type=="surv") {
    group_specific_methods <- c("km", "iptw_km", "iptw_cox",
                                "strat_amato", "strat_nieto")
  } else {
    group_specific_methods <- c("aalen_johansen")
  }

  # keep only events
  event_dat <- mids$data[mids$data[, event]==cause,]

  if (method %in% group_specific_methods) {
    # calculate maximal observed time in each group
    out <- vector(mode="numeric", length=length(levs))

    for (i in seq_len(length(levs))) {
      dat_i <- event_dat[event_dat[, variable]==levs[i],]

      if (nrow(dat_i)==0) {
        max_t <- 0
      } else {
        max_t <- max(dat_i[, ev_time], na.rm=TRUE)
      }

      out[[i]] <- max_t
    }
  } else {
    out <- max(event_dat[, ev_time], na.rm=TRUE)
  }

  max_t_group <- data.frame(group=levs, max_t=out)

  return(max_t_group)
}
