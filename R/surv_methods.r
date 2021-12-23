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

# Assign to global, to get rid off devtools::check() note
utils::globalVariables(c("gaussian", "id"))

## S3 print method for adjustedsurv.method objects
#' @export
print.adjustedsurv.method <- function(x, ...) {
  print(x$plotdata, ...)
}

## S3 summary method for adjustedsurv.method objects
#' @export
summary.adjustedsurv.method <- function(object, ...) {
  summary(object$plotdata, ...)
}

## simple Kaplan-Meier estimate
#' @export
surv_km <- function(data, variable, ev_time, event, conf_int,
                    conf_level=0.95, times=NULL, conf_type="log") {

  form <- paste0("survival::Surv(", ev_time, ", ", event, ") ~ ", variable)

  surv <- survival::survfit.formula(stats::as.formula(form), data=data,
                                    se.fit=conf_int, conf.int=conf_level,
                                    conf.type=conf_type)
  plotdata <- data.frame(time=surv$time,
                         surv=surv$surv)
  # get grouping variable
  group <- c()
  for (strat in names(surv$strata)) {
    group <- c(group, rep(strat, surv$strata[strat]))
  }
  group <- gsub(paste0(variable, "="), "", group)
  plotdata$group <- group

  # get se and confidence interval
  if (conf_int) {
    plotdata$se <- surv$std.err
    plotdata$ci_lower <- surv$lower
    plotdata$ci_upper <- surv$upper
  }

  if (!is.null(times)) {
    plotdata <- specific_times(plotdata, times)
  # ensure that it starts with 0
  } else if (!0 %in% plotdata$time) {
    if (conf_int) {
      levs <- unique(data[, variable])

      row_0 <- data.frame(time=0, group=NA, surv=1, se=0, ci_lower=1,
                          ci_upper=1)
      row_0 <- row_0[rep(1, each=length(levs)), ]
      row_0$group <- levs
      rownames(row_0) <- NULL
      plotdata <- rbind(row_0, plotdata)
    }
  }

  output <- list(plotdata=plotdata,
                 survfit_object=surv)
  class(output) <- "adjustedsurv.method"

  return(output)
}

## IPTW Kaplan-Meier estimate
#' @export
surv_iptw_km <- function(data, variable, ev_time, event, conf_int,
                         conf_level=0.95, times=NULL, treatment_model,
                         weight_method="ps", stabilize=TRUE,
                         trim=FALSE, ...) {

  levs <- levels(data[, variable])

  # get weights
  if (is.numeric(treatment_model)) {
    weights <- treatment_model
    weights <- trim_weights(weights=weights, trim=trim)
  } else {
    weights <- get_iptw_weights(data=data, treatment_model=treatment_model,
                                weight_method=weight_method,
                                variable=variable, stabilize=stabilize,
                                trim=trim, ...)
  }

  plotdata <- vector(mode="list", length=length(levs))
  for (i in seq_len(length(levs))){

    dat_group <- data[data[, variable]==levs[i], ]
    weights_group <- weights[data[, variable]==levs[i]]

    # calculate weighted risk set and events
    tj <- c(0, sort(unique(dat_group[, ev_time][dat_group[, event]==1])))
    dj <- vapply(tj, function(x){sum(weights_group[dat_group[, ev_time]==x &
                                                   dat_group[, event]==1])},
                 FUN.VALUE=numeric(1))
    yj <- vapply(tj, function(x){sum(weights_group[dat_group[, ev_time]>=x])},
                 FUN.VALUE=numeric(1))
    st <- cumprod(1 - (dj / yj))

    plotdata_group <- data.frame(time=tj, surv=st, group=levs[i])

    # adding approximate confidence intervals
    if (conf_int) {

      # variance calculation based of Xie et al.
      m <- vapply(tj, FUN.VALUE=numeric(1),
                  function(x){sum((weights_group[dat_group[, ev_time]>=x])^2)})
      mj <- ((yj^2) / m)
      fracj <- dj / (mj * (yj - dj))
      fracj_cumsum <- cumsum(fracj)
      Vst <- (st^2) * fracj_cumsum

      # standard error
      # NOTE: the n is already accounted for, so taking the sqrt() is
      #       enough to calculate the standard error
      plotdata_group$se <- sqrt(Vst)

      # normal approximation
      surv_ci <- confint_surv(surv=plotdata_group$surv,
                              se=plotdata_group$se,
                              conf_level=conf_level,
                              conf_type="plain")
      plotdata_group$ci_lower <- surv_ci$left
      plotdata_group$ci_upper <- surv_ci$right

    }
    plotdata[[i]] <- plotdata_group
  }

  plotdata <- dplyr::bind_rows(plotdata)

  if (!is.null(times)) {
    plotdata <- specific_times(plotdata, times)
  }

  output <- list(plotdata=plotdata,
                 weights=weights)
  class(output) <- "adjustedsurv.method"

  return(output)
}

## IPTW with univariate cox-model
# NOTE:
# - using the G-Formula directly on the iptw cox model does not work,
#   probably because all predict() methods ignore the weights when
#   calculating the baseline hazard. Only this version is actually unbiased
#' @export
surv_iptw_cox <- function(data, variable, ev_time, event, conf_int,
                          conf_level=0.95, times=NULL, treatment_model,
                          weight_method="ps", stabilize=TRUE,
                          trim=FALSE, ...) {
  # get weights
  if (is.numeric(treatment_model)) {
    weights <- treatment_model
  } else {
    weights <- get_iptw_weights(data=data, treatment_model=treatment_model,
                                weight_method=weight_method,
                                variable=variable, stabilize=stabilize,
                                trim=trim, ...)
  }

  # univariate, weighted cox model
  form <- paste0("survival::Surv(", ev_time, ", ", event, ") ~ strata(",
                 variable, ")")
  model <- survival::coxph(stats::as.formula(form), weights=weights, data=data)
  surv <- survival::survfit(model, se.fit=conf_int, conf.int=conf_level)

  plotdata <- data.frame(time=surv$time,
                         surv=surv$surv)

  # get grouping variable
  group <- c()
  for (strat in names(surv$strata)) {
    group <- c(group, rep(strat, surv$strata[strat]))
  }
  group <- gsub(paste0(variable, "="), "", group)
  plotdata$group <- group

  # get se and confidence interval
  if (conf_int) {
    plotdata$se <- surv$std.err
    plotdata$ci_lower <- surv$lower
    plotdata$ci_upper <- surv$upper
  }

  if (!is.null(times)) {
    plotdata <- specific_times(plotdata, times)
  }

  output <- list(plotdata=plotdata,
                 cox_model=model,
                 survfit_object=surv,
                 weights=weights)
  class(output) <- "adjustedsurv.method"

  return(output)
}

## Using Pseudo-Observations and IPTW
#' @export
surv_iptw_pseudo <- function(data, variable, ev_time, event, conf_int,
                             conf_level=0.95, times, treatment_model,
                             weight_method="ps", stabilize=TRUE,
                             trim=FALSE, se_method="cochrane",
                             censoring_vars=NULL, ipcw_method="binder", ...) {
  levs <- levels(data[, variable])

  # get weights
  if (is.numeric(treatment_model)) {
    weights <- treatment_model
    weights <- trim_weights(weights, trim)
  } else {
    weights <- get_iptw_weights(data=data, treatment_model=treatment_model,
                                weight_method=weight_method,
                                variable=variable, stabilize=stabilize,
                                trim=trim, ...)
  }

  # estimate pseudo observations
  pseudo <- calc_pseudo_surv(data=data,
                             ev_time=ev_time,
                             event=event,
                             times=times,
                             censoring_vars=censoring_vars,
                             ipcw.method=ipcw_method)

  # take weighted mean
  levs <- levels(data[, variable])
  plotdata <- vector(mode="list", length=length(levs))
  for (i in seq_len(length(levs))) {
    surv_lev <- pseudo[data[, variable]==levs[i], ]
    surv_lev <- apply(surv_lev, 2, stats::weighted.mean,
                      w=weights[data[, variable]==levs[i]],
                      na.rm=TRUE)

    data_temp <- data.frame(time=times, surv=surv_lev, group=levs[i])

    if (conf_int) {
      # approximate variance by calculating a
      # weighted version of the variance
      surv_lev <- pseudo[data[, variable]==levs[i], ]

      surv_sd <- apply(surv_lev, 2, weighted.var.se,
                       w=weights[data[, variable]==levs[i]],
                       na.rm=TRUE, se_method=se_method)
      data_temp$se <- sqrt(surv_sd)

      surv_cis <- confint_surv(surv=data_temp$surv, se=data_temp$se,
                               conf_level=conf_level, conf_type="plain")
      data_temp$ci_lower <- surv_cis$left
      data_temp$ci_upper <- surv_cis$right

    }
    plotdata[[i]] <- data_temp
  }
  plotdata <- dplyr::bind_rows(plotdata)
  rownames(plotdata) <- NULL

  output <- list(plotdata=plotdata,
                 pseudo_values=pseudo,
                 weights=weights)
  class(output) <- "adjustedsurv.method"

  return(output)
}

## Direct Adjustment
#' @export
surv_direct <- function(data, variable, ev_time, event, conf_int,
                        conf_level=0.95, times, outcome_model,
                        verbose=FALSE, predict_fun=NULL,  ...) {

  # using the ate function
  if (inherits(outcome_model, c("coxph", "cph")) & is.null(predict_fun)) {

    # NOTE: This prevents an error message in predictSurv() when
    #       the 'weights' argument was used in the coxph() call.
    #       Standard error calculations might be off in that case, not sure.
    outcome_model$weights <- 1
    outcome_model$naive.var <- NULL

    surv <- riskRegression::ate(event=outcome_model, treatment=variable,
                                data=data, estimator="Gformula",
                                times=times, se=conf_int, verbose=verbose,
                                cause=1, ...)
    plotdata <- data.frame(time=surv$meanRisk$time,
                           surv=1 - surv$meanRisk$estimate,
                           group=surv$meanRisk$treatment)

    if (conf_int) {
      plotdata$se <- surv$meanRisk$se

      confint.ate <- utils::getFromNamespace("confint.ate", "riskRegression")

      cis <- confint.ate(surv, level=conf_level)$meanRisk
      plotdata$ci_lower <- 1 - cis$upper
      plotdata$ci_upper <- 1 - cis$lower
    }

    output <- list(plotdata=plotdata,
                   ate_object=surv)
    class(output) <- "adjustedsurv.method"
  # using outcome_model specific functions
  } else {
    plotdata <- surv_g_comp(outcome_model=outcome_model,
                            data=data,
                            variable=variable,
                            times=times,
                            predict_fun=predict_fun,
                            ...)
    output <- list(plotdata=plotdata)
    class(output) <- "adjustedsurv.method"
  }

  return(output)
}

## using models other than coxph in method="direct" for
## survival endpoints
surv_g_comp <- function(outcome_model, data, variable, times,
                        predict_fun, ...) {

  # perform G-Computation
  levs <- levels(data[, variable])
  data_temp <- data
  plotdata <- vector(mode="list", length=length(levs))
  for (i in seq_len(length(levs))) {

    # set variable to one level each
    data_temp[, variable] <- factor(levs[i], levels=levs)

    # use user-supplied custom prediction function
    if (!is.null(predict_fun)) {
      surv_lev <- predict_fun(outcome_model,
                              newdata=data_temp,
                              times=times,
                              ...)
    # use function from pec package for fast & easy survival prediction
    # NOTE: while aalen & cox.aalen are not working in predictRisk, keep them
    #       here
    } else if (inherits(outcome_model, c("pecCforest", "pecRpart",
                                         "selectCox", "aalen", "cox.aalen"))) {
      requireNamespace("pec")

      # get model.matrix if needed
      if (inherits(outcome_model, c("aalen"))) {
        mod_vars <- all.vars(outcome_model$call$formula)
        mod_form <- paste0(" ~ ", paste0(mod_vars, collapse=" + "))
        mod_data <- as.data.frame(stats::model.matrix(
          stats::as.formula(mod_form), data=data_temp))
      } else {
        mod_data <- data_temp
      }

      # predict survival
      surv_lev <- pec::predictSurvProb(outcome_model,
                                       newdata=mod_data,
                                       times=times,
                                       ...)
    # using predictRisk
    # NOTE: - In this context "BinaryTree", "lrm", "rpart" make no sense
    #         because they don't allow prediction at multiple t
    #       - Can't use "coxph.penal" due to bugs in predictRisk
    #       - "hal9001" not tested, "singleEventCB" has problems
    } else if (inherits(outcome_model, c("coxphTD", "prodlim", "psm",
                                         "ranger", "rfsrc", "riskRegression",
                                         "ARR", "penfitS3", "gbm",
                                         "flexsurvreg", "singleEventCB",
                                         "wglm", "hal9001"))) {
      requireNamespace("riskRegression")

      surv_lev <- riskRegression::predictRisk(object=outcome_model,
                                              newdata=data_temp,
                                              times=times,
                                              cause=1,
                                              ...)
      surv_lev <- 1 - surv_lev

    # using predictProb
    } else if (inherits(outcome_model, c("glm", "ols", "randomForest"))) {
      requireNamespace("pec")
      predictProb <- utils::getFromNamespace("predictProb", "pec")

      surv_lev <- predictProb(object=outcome_model,
                              newdata=data_temp,
                              times=times,
                              ...)
    # some customly created functions from here
    } else if (inherits(outcome_model, "mexhaz")) {
      # can't do predictions at 0
      times <- times[times > 0]
      times <- c(0.000001, times)
      surv_lev <- vapply(X=times, object=outcome_model, data.val=data_temp,
                   FUN=function(x, ...){
                     stats::predict(time.pts=x, ...)$results$surv},
                   FUN.VALUE=numeric(nrow(data)))
    # try to directly use S3 prediction function
    } else {
      surv_lev <- tryCatch(
        expr={stats::predict(outcome_model,
                             newdata=data_temp,
                             times=times,
                             ...)},
        error=function(e){stop("The following error occured using",
                               " the default S3 predict method: '", e,
                               "' Specify a valid 'predict_fun' or",
                               " use a different model. See details.")}
      )
    }

    # take arithmetic mean of predictions and add those to the
    # output object
    surv_lev <- apply(X=surv_lev, MARGIN=2, FUN=mean, na.rm=TRUE)
    row <- data.frame(time=times, surv=surv_lev, group=levs[i])
    plotdata[[i]] <- row

  }
  plotdata <- dplyr::bind_rows(plotdata)
  row.names(plotdata) <- seq_len(nrow(plotdata))

  return(plotdata)
}

## Using Propensity Score Matching
#' @export
surv_matching <- function(data, variable, ev_time, event, conf_int=FALSE,
                          conf_level=0.95, times, treatment_model,
                          gtol=0.001, ...) {

  # if it's a factor, turn it into numeric
  if (is.factor(data[, variable])) {
    levs <- levels(data[, variable])
    data[, variable] <- ifelse(data[,variable]==levs[1], 0, 1)
  } else {
    levs <- sort(unique(data[, variable]))
  }

  if (is.numeric(treatment_model)) {
    ps_score <- treatment_model
  } else {
    ps_score <- stats::predict.glm(treatment_model, newdata=data,
                                   type="response")
  }

  # trim extreme propensity score according to gtol
  ps_score[ps_score < gtol] <- gtol
  ps_score[ps_score > (1 - gtol)] <- 1 - gtol

  # perform matching
  rr <- Matching::Match(Tr=data[, variable], X=ps_score, estimand="ATE", ...)
  m_dat <- rbind(data[rr$index.treated, ], data[rr$index.control, ])

  weights <- rr$weights
  m_dat$match_weights <- c(weights, weights)

  # estimate survival curve
  form <- paste0("survival::Surv(", ev_time, ", ", event, ") ~ ", variable)
  surv <- survival::survfit(stats::as.formula(form), data=m_dat,
                            se.fit=conf_int, conf.int=conf_level,
                            weights=m_dat$match_weights,
                            robust=TRUE)
  plotdata <- data.frame(time=surv$time,
                         surv=surv$surv,
                         group=c(rep(levs[1], surv$strata[1]),
                                 rep(levs[2], surv$strata[2])))

  if (!is.null(times)) {
    plotdata <- specific_times(plotdata, times)
  }

  output <- list(plotdata=plotdata,
                 match_object=rr,
                 survfit_object=surv)
  class(output) <- "adjustedsurv.method"

  return(output)
}

## Using Augmented Inverse Probability of Treatment Weighting
#' @export
surv_aiptw <- function(data, variable, ev_time, event, conf_int,
                       conf_level=0.95, times, outcome_model,
                       treatment_model, censoring_model=NULL,
                       verbose=FALSE, ...) {

  # defaults for input models
  if (is.null(censoring_model)) {
    form <- paste0("survival::Surv(", ev_time, ", ", event, "==0) ~ 1")
    censoring_model <- survival::coxph(stats::as.formula(form), data=data,
                                       x=TRUE)
  }

  # estimate AIPTW cumulative incidence
  curve <- riskRegression::ate(event=outcome_model,
                               treatment=treatment_model,
                               censor=censoring_model,
                               data=data,
                               times=times,
                               se=conf_int,
                               verbose=verbose,
                               estimator="AIPTW,AIPCW",
                               cause=1)
  # calculate survival estimates
  plotdata <- data.frame(time=curve$meanRisk$time,
                         surv=1 - curve$meanRisk$estimate,
                         group=curve$meanRisk$treatment)
  if (conf_int) {
    plotdata$se <- curve$meanRisk$se

    confint.ate <- utils::getFromNamespace("confint.ate", "riskRegression")

    cis <- confint.ate(curve, level=conf_level, ci=TRUE)$meanRisk
    plotdata$ci_lower <- 1 - cis$upper
    plotdata$ci_upper <- 1 - cis$lower
  }

  output <- list(plotdata=plotdata,
                 ate_object=curve)
  class(output) <- "adjustedsurv.method"

  return(output)
}

## Using Pseudo Observations and Direct Adjustment
#' @export
surv_direct_pseudo <- function(data, variable, ev_time, event,
                               conf_int=FALSE, conf_level=0.95, times,
                               outcome_vars, type_time="factor",
                               spline_df=5, censoring_vars=NULL,
                               ipcw_method="binder") {

  # estimate pseudo observations
  pseudo <- calc_pseudo_surv(data=data,
                             ev_time=ev_time,
                             event=event,
                             times=times,
                             censoring_vars=censoring_vars,
                             ipcw.method=ipcw_method)

  # remove "variable" from outcome_vars because it is always included
  outcome_vars <- outcome_vars[outcome_vars!=variable]

  # some constants
  len <- length(times)
  n <- nrow(data)
  group <- data[, variable]

  # create data for geese
  Sdata <- data.frame(yi=1 - c(pseudo),
                      group=rep(group, len),
                      vtime=rep(times, rep(n, len)),
                      id=rep(1:n, len))
  for (col in outcome_vars) {
    Sdata[, col] <- rep(data[, col], len)
  }

  if (type_time=="factor") {
    Sdata$vtime <- as.factor(Sdata$vtime)

    if (length(times)==1) {
      geese_formula <- paste("yi ~ ", paste(outcome_vars, collapse=" + "),
                             " + group")
    } else {
      geese_formula <- paste("yi ~ vtime + ",
                             paste(outcome_vars, collapse=" + "),
                             " + group")
    }

  } else if (type_time=="bs") {
    geese_formula <- paste("yi ~ splines::bs(vtime, df=", spline_df, ") + ",
                           paste(outcome_vars, collapse=" + "), " + group")
  } else if (type_time=="ns") {
    geese_formula <- paste("yi ~ splines::ns(vtime, df=", spline_df, ") + ",
                           paste(outcome_vars, collapse=" + "), " + group")
  }

  # remove rows where pseudo-values are NA for geese
  Sdata_fit <- Sdata[!is.na(Sdata$yi), ]

  # call geese
  geese_mod <- geepack::geese(stats::as.formula(geese_formula), scale.fix=TRUE,
                              data=Sdata_fit, family=gaussian, id=id,
                              jack=FALSE, mean.link="cloglog",
                              corstr="independence")

  # initialize outcome df list
  levs <- levels(data[, variable])
  plotdata <- vector(mode="list", length=length(levs))

  # do direct adjustment
  for (i in seq_len(length(levs))) {

    Sdata$group <- factor(levs[i], levels=levs)
    pred <- geese_predictions(geese_mod, Sdata, times=times, n=n)

    m <- exp(-exp(pred))
    surv <- apply(m, 2, mean, na.rm=TRUE)

    plotdata[[i]] <- data.frame(time=times, surv=surv, group=levs[i])

  }
  plotdata <- dplyr::bind_rows(plotdata)
  rownames(plotdata) <- NULL

  output <- list(plotdata=plotdata,
                 pseudo_values=pseudo,
                 geese_model=geese_mod)
  class(output) <- "adjustedsurv.method"

  return(output)
}

## Using Pseudo Observations Doubly-Robust, Wang (2018)
#' @export
surv_aiptw_pseudo <- function(data, variable, ev_time, event, conf_int,
                              conf_level=0.95, times, outcome_vars,
                              treatment_model, type_time="factor",
                              spline_df=5, censoring_vars=NULL,
                              ipcw_method="binder") {
  # some constants
  len <- length(times)
  n <- nrow(data)
  group <- data[, variable]

  if (is.numeric(treatment_model)) {
    ps_score <- treatment_model
  } else if (inherits(treatment_model, "glm")) {
    ps_score <- treatment_model$fitted.values
  } else if (inherits(treatment_model, "multinom")) {
    predict.multinom <- utils::getFromNamespace("predict.multinom", "nnet")
    ps_score <- predict.multinom(treatment_model, newdata=data, type="probs")
  }

  # estimate pseudo observations
  pseudo <- calc_pseudo_surv(data=data,
                             ev_time=ev_time,
                             event=event,
                             times=times,
                             censoring_vars=censoring_vars,
                             ipcw.method=ipcw_method)
  # create data for geese
  Sdata <- data.frame(yi=c(pseudo),
                      group=rep(group, len),
                      vtime=rep(times, rep(n, len)),
                      id=rep(1:n, len))

  outcome_vars <- outcome_vars[outcome_vars != variable]
  for (col in outcome_vars) {
    Sdata[, col] <- rep(data[, col], len)
  }

  if (type_time=="factor") {
    Sdata$vtime <- as.factor(Sdata$vtime)
    geese_formula <- paste("yi ~ vtime + ", paste(outcome_vars, collapse=" + "),
                           " + group")
  } else if (type_time=="bs") {
    geese_formula <- paste("yi ~ splines::bs(vtime, df=", spline_df, ") + ",
                           paste(outcome_vars, collapse=" + "), " + group")
  } else if (type_time=="ns") {
    geese_formula <- paste("yi ~ splines::ns(vtime, df=", spline_df, ") + ",
                           paste(outcome_vars, collapse=" + "), " + group")
  }

  # remove rows where pseudo-values are NA for geese
  Sdata_fit <- Sdata[!is.na(Sdata$yi), ]

  # call geese
  geese_mod <- geepack::geese(stats::as.formula(geese_formula), scale.fix=TRUE,
                              data=Sdata_fit, family=gaussian, id=id,
                              jack=FALSE, mean.link="cloglog",
                              corstr="independence")

  # get direct adjustment estimates
  levs <- levels(data[, variable])
  plotdata <- vector(mode="list", length=length(levs))

  for (i in seq_len(length(levs))) {

    Sdata$group <- factor(levs[i], levels=levs)
    pred <- geese_predictions(geese_mod, Sdata, times=times, n=n)

    m <- 1 - exp(-exp(pred))

    # augment estimates using propensity score
    if (length(levs) > 2) {
      group_ind <- ifelse(group==levs[i], 1, 0)
      ps_score_lev <- ps_score[, levs[i]]

      dr <- (pseudo*group_ind-(group_ind-ps_score_lev)*m)/ps_score_lev
    # if binary, use equation from the paper directly
    } else if (i == 1) {
      group <- ifelse(data[, variable]==levs[1], 0, 1)
      dr <- (pseudo*(1-group)+(group-ps_score)*m)/(1-ps_score)
    } else if (i == 2) {
      group <- ifelse(data[, variable]==levs[1], 0, 1)
      dr <- (pseudo*group-(group-ps_score)*m)/ps_score
    }

    surv <- apply(dr, 2, mean, na.rm=TRUE)

    if (conf_int) {

      pseudo_dr_se <- function(x, n, na.rm) {
        sqrt(stats::var(x, na.rm=na.rm) / n)
      }
      surv_se <- apply(dr, 2, pseudo_dr_se, n=n, na.rm=TRUE)

      surv_ci <- confint_surv(surv=surv, se=surv_se, conf_level=conf_level,
                              conf_type="plain")

      plotdata[[i]] <- data.frame(time=times, surv=surv, group=levs[i],
                                  se=surv_se, ci_lower=surv_ci$left,
                                  ci_upper=surv_ci$right)

    } else {
      plotdata[[i]] <- data.frame(time=times, surv=surv, group=levs[i])
    }

  }
  plotdata <- dplyr::bind_rows(plotdata)
  rownames(plotdata) <- NULL

  output <- list(plotdata=plotdata,
                 pseudo_values=pseudo,
                 geese_model=geese_mod)
  class(output) <- "adjustedsurv.method"

  return(output)
}

## Using Empirical Likelihood Estimation
#' @export
surv_emp_lik <- function(data, variable, ev_time, event, conf_int=FALSE,
                         times, treatment_vars, moment="first",
                         standardize=FALSE, gtol=0.00001,
                         max_iter=100, newton_tol=1.0e-06) {

  # if it's a factor, turn it into numeric
  if (is.factor(data[, variable])) {
    levs <- levels(data[, variable])
    data[, variable] <- ifelse(data[, variable]==levs[1], 0, 1)
  } else {
    levs <- sort(unique(data[, variable]))
  }

  el_0 <- el.est(y=data[, ev_time],
                 delta=data[, event],
                 treat=data[, variable],
                 x=as.matrix(data[, treatment_vars]),
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
                 x=as.matrix(data[, treatment_vars]),
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

  output <- list(plotdata=plotdata)
  class(output) <- "adjustedsurv.method"

  return(output)
}

## Targeted Maximum Likelihood Estimation
#' @export
surv_tmle <- function(data, variable, ev_time, event, conf_int,
                      conf_level=0.95, times, adjust_vars=NULL,
                      SL.ftime=NULL, SL.ctime=NULL, SL.trt=NULL,
                      glm.ftime=NULL, glm.ctime=NULL, glm.trt=NULL,
                      ...) {

  # if it's a factor, turn it into numeric
  if (is.factor(data[, variable])) {
    levs <- levels(data[, variable])
    data[,variable] <- ifelse(data[, variable]==levs[1], 0, 1)
  } else {
    levs <- unique(data[, variable])
  }

  # gather needed data
  if (is.null(adjust_vars)) {
    all_covars <- colnames(data)
    all_covars <- all_covars[!all_covars %in% c(variable, ev_time, event)]
    adjust_vars <- data[, all_covars]
  } else {
    adjust_vars <- data[, adjust_vars]
  }

  # TMLE fit
  fit <- survtmle::survtmle(
    ftime=data[, ev_time],
    ftype=data[, event],
    trt=data[, variable],
    t0=max(data[, ev_time]),
    adjustVars=adjust_vars,
    method="hazard",
    returnIC=TRUE,
    verbose=FALSE,
    returnModels=TRUE,
    SL.trt=SL.trt,
    SL.ftime=SL.ftime,
    SL.ctime=SL.ctime,
    glm.ftime=glm.ftime,
    glm.ctime=glm.ctime,
    glm.trt=glm.trt,
    ...
  )
  # extract cumulative incidence at each timepoint overall
  tpfit <- withCallingHandlers({
    survtmle.timepoints(fit, times=times, SL.trt=SL.trt, SL.ctime=SL.ctime,
                        SL.ftime=SL.ftime, glm.trt=glm.trt,
                        glm.ctime=glm.ctime, glm.ftime=glm.ftime)
  }, warning=function(w) {
    if (startsWith(conditionMessage(w), "Using formula(x) is deprecated"))
      invokeRestart("muffleWarning")
  })

  s_0 <- 1 - unlist(lapply(tpfit, function(x) {x$est[1]}))
  s_1 <- 1 - unlist(lapply(tpfit, function(x) {x$est[2]}))

  # put together
  plotdata <- data.frame(time=rep(times, 2),
                         surv=c(s_0, s_1),
                         group=c(rep(levs[1], length(times)),
                                 rep(levs[2], length(times))))

  if (conf_int) {

    var_0 <- unlist(lapply(tpfit, function(x) {x$var[1]}))
    var_1 <- unlist(lapply(tpfit, function(x) {x$var[2]}))

    plotdata$se <- sqrt(c(var_0, var_1))

    confint.tp.survtmle <- utils::getFromNamespace("confint.tp.survtmle",
                                                   "survtmle")
    survtmle_ci <- confint.tp.survtmle(tpfit, level=conf_level)

    plotdata$ci_lower <- 1 - c(survtmle_ci$`0 1`[, 1], survtmle_ci$`1 1`[, 1])
    plotdata$ci_upper <- 1 - c(survtmle_ci$`0 1`[, 2], survtmle_ci$`1 1`[, 2])

  }

  output <- list(plotdata=plotdata,
                 survtmle_object=fit,
                 survtmle.timepoints_object=tpfit)
  class(output) <- "adjustedsurv.method"

  return(output)
}

## One-Step Targeted Maximum Likelihood Estimation
#' @export
surv_ostmle <- function(data, variable, ev_time, event, conf_int,
                        conf_level=0.95, times, adjust_vars=NULL,
                        SL.ftime=NULL, SL.ctime=NULL, SL.trt=NULL,
                        epsilon=1, max_num_iteration=100,
                        psi_moss_method="l2", tmle_tolerance=NULL,
                        gtol=1e-3) {

  # if it's a factor, turn it into numeric
  if (is.factor(data[, variable])) {
    levs <- levels(data[, variable])
    data[, variable] <- ifelse(data[, variable]==levs[1], 0, 1)
  } else {
    levs <- unique(data[, variable])
  }

  # gather needed data
  if (is.null(adjust_vars)) {
    all_covars <- colnames(data)
    all_covars <- all_covars[!all_covars %in% c(variable, ev_time, event)]
    adjust_vars <- data[, all_covars]
  } else {
    adjust_vars <- data[, adjust_vars]
  }

  # time point grid
  k_grid <- 1:max(data[, ev_time])

  # create initial fit object
  sl_fit <- initial_sl_fit(
    T_tilde=data[, ev_time],
    Delta=data[, event],
    A=data[, variable],
    W=adjust_vars,
    t_max=max(data[, ev_time]),
    sl_treatment=SL.trt,
    sl_censoring=SL.ctime,
    sl_failure=SL.ftime,
    gtol=gtol
  )

  # set to same time point grid
  sl_fit$density_failure_1$t <- k_grid
  sl_fit$density_failure_0$t <- k_grid
  sl_fit$density_censor_1$t <- k_grid
  sl_fit$density_censor_0$t <- k_grid

  invisible(sl_fit$density_failure_1$hazard_to_survival())
  invisible(sl_fit$density_failure_0$hazard_to_survival())

  # for treatment == 1
  moss_hazard_fit_1 <- MOSS_hazard$new(
    A=data[, variable],
    T_tilde=data[, ev_time],
    Delta=data[, event],
    density_failure=sl_fit$density_failure_1,
    density_censor=sl_fit$density_censor_1,
    g1W=sl_fit$g1W,
    A_intervene=1,
    k_grid=k_grid
  )
  psi_moss_hazard_1 <- moss_hazard_fit_1$iterate_onestep(
    epsilon=epsilon,
    max_num_interation=max_num_iteration,
    tmle_tolerance=tmle_tolerance,
    verbose=FALSE,
    method=psi_moss_method
  )

  # for treatment == 0
  moss_hazard_fit_0 <- MOSS_hazard$new(
    A=data[, variable],
    T_tilde=data[, ev_time],
    Delta=data[, event],
    density_failure=sl_fit$density_failure_0,
    density_censor=sl_fit$density_censor_0,
    g1W=sl_fit$g1W,
    A_intervene=0,
    k_grid=k_grid
  )
  psi_moss_hazard_0 <- moss_hazard_fit_0$iterate_onestep(
    epsilon=epsilon,
    max_num_interation=max_num_iteration,
    tmle_tolerance=tmle_tolerance,
    verbose=FALSE,
    method=psi_moss_method
  )

  # put together
  plotdata <- data.frame(time=c(k_grid, k_grid),
                         surv=c(psi_moss_hazard_0, psi_moss_hazard_1),
                         group=c(rep(levs[1], length(k_grid)),
                                 rep(levs[2], length(k_grid))))
  # keep only time points in times
  plotdata <- plotdata[which(plotdata$time %in% times), ]

  if (conf_int) {

    # for treatment == 1
    s_1 <- as.vector(psi_moss_hazard_1)
    eic_fit <- eic$new(
      A = data[, variable],
      T_tilde = data[, ev_time],
      Delta = data[, event],
      density_failure = moss_hazard_fit_1$density_failure,
      density_censor = moss_hazard_fit_1$density_censor,
      g1W = moss_hazard_fit_1$g1W,
      psi = s_1,
      A_intervene = 1
    )
    eic_matrix <- eic_fit$all_t(k_grid=k_grid)
    se_1 <- compute_se_moss(eic_matrix, alpha=1-conf_level)


    # for treatment == 0
    s_0 <- as.vector(psi_moss_hazard_0)
    eic_fit <- eic$new(
      A = data[, variable],
      T_tilde = data[, ev_time],
      Delta = data[, event],
      density_failure = moss_hazard_fit_0$density_failure,
      density_censor = moss_hazard_fit_0$density_censor,
      g1W = moss_hazard_fit_0$g1W,
      psi = s_0,
      A_intervene = 1
    )
    eic_matrix <- eic_fit$all_t(k_grid=k_grid)
    se_0 <- compute_se_moss(eic_matrix, alpha=1-conf_level)

    plotdata$se <- c(se_0, se_1)

    crit <- stats::qnorm(1-((1-conf_level)/2))
    plotdata$ci_lower <- plotdata$surv - crit * plotdata$se
    plotdata$ci_upper <- plotdata$surv + crit * plotdata$se

  }

  output <- list(plotdata=plotdata,
                 psi_moss_hazard_0=psi_moss_hazard_0,
                 psi_moss_hazard_1=psi_moss_hazard_1)
  class(output) <- "adjustedsurv.method"

  return(output)
}

## Adjustment based on a weighted average of stratified Kaplan-Meier estimates
## using the method by Cupples et al.
#' @export
surv_strat_cupples <- function(data, variable, ev_time, event,
                               conf_int=FALSE, conf_level=0.95, times,
                               adjust_vars, reference=NULL, na.rm=FALSE) {

  # devtools checks
  . <- .COVARS <- time <- treat_group <- surv <- count <- NULL

  # additional variables needed
  # NOTE: This code assumes that there is no column named .ALL or .COVARS
  #       and that there are no tabs in the column names
  data$.ALL <- interaction(data[, c(variable, adjust_vars)], sep="\t")
  if (is.null(reference)) {
    reference <- data
  }
  reference$.COVARS <- interaction(reference[, adjust_vars], sep="\t")

  # Kaplan-Meier survival curve for each possible strata at
  # every event time
  plotdata <- surv_km(data=data,
                      variable=".ALL",
                      ev_time=ev_time,
                      event=event,
                      times=times,
                      conf_int=FALSE)$plotdata

  # add indicator for treatment group and remove said treatment group from
  # the 'group' variable
  plotdata$treat_group <- sub("\t.*", "", plotdata$group)
  plotdata$group <- sub(".*?\t", "", plotdata$group)

  # add count of 'group' at baseline, ignoring treatment group
  count_dat <- reference %>%
    dplyr::group_by(.COVARS) %>%
    dplyr::summarise(count = dplyr::n())
  colnames(count_dat)[1] <- c("group")

  # merge together
  plotdata <- merge(count_dat, plotdata, by="group", all.y=TRUE)

  # take weighted in each treatment group at each t, over strata
  plotdata <- plotdata %>%
    dplyr::group_by(., time, treat_group) %>%
    dplyr::summarise(surv=stats::weighted.mean(x=surv, w=count, na.rm=FALSE),
                     .groups="drop_last")
  colnames(plotdata) <- c("time", "group", "surv")

  # remove NAs
  if (na.rm) {
    plotdata <- plotdata[!is.na(plotdata$surv), ]
  }

  output <- list(plotdata=as.data.frame(plotdata))
  class(output) <- "adjustedsurv.method"

  return(output)
}

## Adjustment based on a weighted average of stratified Kaplan-Meier estimates
## using the method by Amato (1988)
#' @export
surv_strat_amato <- function(data, variable, ev_time, event,
                             conf_int=FALSE, conf_level=0.95,
                             times=NULL, adjust_vars, reference=NULL) {

  # silence checks
  . <- group <- time <- wdj <- wcj <- wrj <- delta_wdj <- delta_wd <- wr <- NULL
  times_input <- times

  # proportions in reference data
  if (is.null(reference)) {
    reference <- data
  }
  reference$.COVARS <- interaction(reference[, adjust_vars])
  Pjs <- prop.table(table(reference$.COVARS))

  # also calculate strata variable in data
  data$.COVARS <- interaction(data[, adjust_vars])

  levs <- levels(data[, variable])
  levs_adjust_var <- levels(data$.COVARS)
  out <- list()
  for (i in seq_len(length(levs))) {

    # data for treatment i
    dat_I <- data[data[, variable]==levs[i], ]
    times <- c(0, sort(unique(dat_I[, ev_time][dat_I[, event]==1])))

    for (j in seq_len(length(levs_adjust_var))) {

      # data for treatment i and only strata j
      dat_IJ <- dat_I[dat_I$.COVARS==levs_adjust_var[j], ]

      # weights for these individuals
      ajn <- nrow(dat_I) * Pjs[levs_adjust_var[j]] / nrow(dat_IJ)

      # people observed to fail
      Ndj <- vapply(times, FUN=function(x) {sum(dat_IJ[, ev_time] <= x &
                                                dat_IJ[, event]==1)},
                    FUN.VALUE=numeric(1))
      delta_Ndj <- vapply(times, FUN=function(x) {sum(dat_IJ[, ev_time] == x &
                                                      dat_IJ[, event]==1)},
                          FUN.VALUE=numeric(1))

      # people censored
      Ncj <- vapply(times, FUN=function(x) {sum(dat_IJ[, ev_time] <= x &
                                                dat_IJ[, event]==0)},
                    FUN.VALUE=numeric(1))
      # people at risk
      Nrj <- vapply(times, FUN=function(x) {sum(dat_IJ[, ev_time] >= x)},
                    FUN.VALUE=numeric(1))

      # put together
      temp <- data.frame(time=times, ajn=ajn[[1]],
                         Ndj=Ndj, Ncj=Ncj, Nrj=Nrj, delta_Ndj=delta_Ndj,
                         strata=levs_adjust_var[j], group=levs[i])
      out[[length(out)+1]] <- temp
    }
  }
  dat_stats <- dplyr::bind_rows(out)

  # calculate sums of weights
  dat_stats$wdj <- dat_stats$Ndj * dat_stats$ajn
  dat_stats$wcj <- dat_stats$Ncj * dat_stats$ajn
  dat_stats$wrj <- dat_stats$Nrj * dat_stats$ajn
  dat_stats$delta_wdj <- dat_stats$delta_Ndj * dat_stats$ajn

  # calculate survival probability
  plotdata <- dat_stats %>%
    dplyr::group_by(., group, time) %>%
    dplyr::summarise(wd=sum(wdj),
                     wc=sum(wcj),
                     wr=sum(wrj),
                     delta_wd=sum(delta_wdj),
                     .groups="drop_last") %>%
    dplyr::mutate(., surv=cumprod(1 - (delta_wd / wr)))

  # remove unnecessary variables
  plotdata$wd <- NULL
  plotdata$wc <- NULL
  plotdata$wr <- NULL
  plotdata$delta_wd <- NULL
  plotdata$group <- factor(plotdata$group, levels=levs)
  plotdata <- as.data.frame(plotdata)

  if (!is.null(times_input)) {
    plotdata <- specific_times(plotdata, times_input)
  }

  output <- list(plotdata=plotdata,
                 Pjs=Pjs)
  class(output) <- "adjustedsurv.method"
  return(output)
}

## Adjustment based on a weighted average of stratified Kaplan-Meier estimates
## using the method by Gregory (1988) and Nieto & Coresh (1996)
# NOTE: Equations are due to Nieto & Coresh (1996) because while both
#       methods produce the same results when using the full data as reference,
#       only Nieto's formulation allows the calculation of confidence intervals.
#' @export
surv_strat_gregory_nieto <- function(data, variable, ev_time, event,
                                     conf_int, conf_level=0.95,
                                     times=NULL, adjust_vars, na.rm=FALSE) {

  # silence checks
  . <- time <- group <- frac <- est_var <- wji <- var_j <- NULL

  data$.COVARS <- interaction(data[, adjust_vars])
  times_input <- times

  # needed levels
  levs <- levels(data[, variable])
  levs_adjust_var <- levels(data$.COVARS)

  out <- list()
  for (i in seq_len(length(levs))) {

    dat_X <- data[data[, variable]==levs[i], ]

    # 1.)
    tj <- c(0, sort(unique(dat_X[, ev_time][dat_X[, event]==1])))

    for (j in seq_len(length(levs_adjust_var))) {

      # data at X, Z
      dat_XZ <- data[data[, variable]==levs[i] &
                       data$.COVARS==levs_adjust_var[j], ]

      # 2.)
      nxzj <- vapply(tj, function(x) {sum(dat_XZ[, ev_time]>=x)},
                     FUN.VALUE=numeric(1))
      axzj <- vapply(tj, function(x) {sum(dat_XZ[, ev_time]==x &
                                           dat_XZ[, event]==1)},
                     FUN.VALUE=numeric(1))
      # 3.)
      qxzj <- axzj / nxzj

      # 4.) but modified, calculating n at risk in strata overall instead
      #     of using the control group
      dat_Z <- data[data$.COVARS==levs_adjust_var[j], ]
      nz <- vapply(tj, function(x) {sum(dat_Z[, ev_time]>=x)},
                   FUN.VALUE=numeric(1))

      out[[length(out)+1]] <- data.frame(time=tj, nxzj=nxzj, axzj=axzj,
                                         qxzj=qxzj, nz=nz, group=levs[i],
                                         strata=levs_adjust_var[j])
    }
  }
  out <- dplyr::bind_rows(out)

  # appendix 1
  if (conf_int) {
    # calculate total nj (n at risk) in full data
    tj_overall <- c(0, sort(unique(data[, ev_time][data[, event]==1])))
    nj <- vapply(tj_overall, function(x) {sum(data[, ev_time]>=x)},
                 FUN.VALUE=numeric(1))
    dat_nj <- data.frame(time=tj_overall, nj=nj)

    # merge to previous output
    out <- merge(out, dat_nj, by="time", all.x=TRUE)

    # calculate wji
    out$wji <- out$nz / out$nj

    # calculate the sum needed at the end of equation 7
    out <- out %>%
      dplyr::group_by(., time, group) %>%
      dplyr::mutate(sum_wq=sum(wji * qxzj))

    # stratum specific variance
    out$var_j <- out$wji^2 * (((1 - out$qxzj)*out$qxzj) / out$nxzj) *
      (1 / (1 - out$sum_wq)^2)
  } else {
    out$var_j <- 0
  }

  # 5.) + 6.) and variance from appendix
  plotdata <- out %>%
    dplyr::group_by(., time, group) %>%
    dplyr::summarise(frac= (sum(qxzj * nz)) / sum(nz),
                     est_var=sum(var_j),
                     .groups="drop_last") %>%
    dplyr::group_by(., group) %>%
    dplyr::mutate(surv=cumprod(1 - frac),
                  se=sqrt(cumsum(est_var)))
  plotdata$frac <- NULL
  plotdata <- as.data.frame(plotdata)

  # equation 8 appendix
  if (conf_int) {
    surv_ci <- confint_surv(surv=log(plotdata$surv),
                            se=plotdata$se,
                            conf_level=conf_level,
                            conf_type="plain")
    plotdata$ci_lower <- exp(surv_ci$left)
    plotdata$ci_upper <- exp(surv_ci$right)
  } else {
    plotdata$se <- NULL
    plotdata$est_var <- NULL
  }

  if (!is.null(times_input)) {
    plotdata <- specific_times(plotdata, times_input)
  }

  # remove NAs
  if (na.rm) {
    plotdata <- plotdata[!is.na(plotdata$surv), ]
  }

  output <- list(plotdata=plotdata)
  class(output) <- "adjustedsurv.method"

  return(output)
}
