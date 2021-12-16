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

# A tradeoff had to be made here. As can clearly be seen, some cif_method
# functions are nearly identical to the surv_method equivalents. Instead
# of writing one more general function, we choose to write separate functions
# to make the documentation and usage easier, while making the code itself
# a little messier with unnecessary repetition.

## S3 print method for adjustedcif.method objects
#' @export
print.adjustedcif.method <- function(x, ...) {
  print(x$plotdata, ...)
}

## S3 summary method for adjustedcif.method objects
#' @export
summary.adjustedcif.method <- function(object, ...) {
  summary(object$plotdata, ...)
}

## Aalen-Johansen estimator
#' @export
cif_aalen_johansen <- function(data, variable, ev_time, event, cause,
                               conf_int, conf_level=0.95, times=NULL, ...) {

  cif <- cmprsk::cuminc(ftime=data[, ev_time],
                        fstatus=data[, event],
                        group=data[, variable],
                        ...)

  levs <- unique(data[, variable])
  cif_names <- paste(levs, cause)
  plotdata <- vector(mode="list", length=length(cif_names))
  for (i in seq_len(length(cif_names))) {
    plotdata[[i]] <- data.frame(time=cif[[cif_names[i]]]$time,
                                cif=cif[[cif_names[i]]]$est,
                                group=levs[i],
                                se=cif[[cif_names[i]]]$var)
  }
  plotdata <- as.data.frame(dplyr::bind_rows(plotdata))

  if (conf_int) {
    cif_cis <- confint_surv(surv=plotdata$cif, se=sqrt(plotdata$se),
                            conf_level=conf_level, conf_type="plain")
    plotdata$ci_lower <- cif_cis$left
    plotdata$ci_upper <- cif_cis$right
  } else {
    plotdata$se <- NULL
  }

  if (!is.null(times)) {
    plotdata <- specific_times(plotdata, times, cif=TRUE)
  }

  # remove weird structure in cmprsk::cuminc call
  ids <- seq(2, nrow(plotdata), 2)
  ids <- ids[1:(length(ids-1))]
  plotdata <- plotdata[ids, ]

  output <- list(plotdata=plotdata,
                 cuminc_object=cif)
  class(output) <- "adjustedcif.method"

  return(output)
}

## IPTW
#' @export
cif_iptw <- function(data, variable, ev_time, event, cause, conf_int,
                     conf_level=0.95, times, treatment_model,
                     censoring_model=NULL, verbose=FALSE, ...) {
  # empty censoring model if not specified
  if (is.null(censoring_model)) {
    form <- paste0("survival::Surv(", ev_time, ", ", event, "==0) ~ 1")
    censoring_model <- survival::coxph(stats::as.formula(form), data=data,
                                       x=TRUE, y=TRUE)
  }

  cif <- riskRegression::ate(event=c(ev_time, event),
                             treatment=treatment_model,
                             data=data,
                             estimator="IPTW",
                             times=times,
                             se=conf_int,
                             verbose=verbose,
                             cause=cause,
                             censor=censoring_model,
                             ...)
  plotdata <- data.frame(time=cif$meanRisk$time,
                         cif=cif$meanRisk$estimate,
                         group=cif$meanRisk$treatment)

  if (conf_int) {
    plotdata$se <- cif$meanRisk$se

    confint.ate <- utils::getFromNamespace("confint.ate", "riskRegression")

    cis <- confint.ate(cif, level=conf_level)$meanRisk
    plotdata$ci_lower <- cis$lower
    plotdata$ci_upper <- cis$upper
  }

  output <- list(plotdata=plotdata,
                 ate_object=cif)
  class(output) <- "adjustedcif.method"

  return(output)
}

# IPTW pseudo
#' @export
cif_iptw_pseudo <- function(data, variable, ev_time, event, cause,
                            conf_int, conf_level=0.95, times,
                            treatment_model, weight_method="ps",
                            stabilize=TRUE, trim=FALSE,
                            se_method="cochrane", ...) {
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

  # estimate pseudo observations
  hist_formula <- stats::as.formula(paste("prodlim::Hist(", ev_time, ", ",
                                          event, ") ~ 1"))
  pseudo <- prodlim::jackknife(prodlim::prodlim(hist_formula, data=data),
                               times=times, cause=cause)

  # take weighted mean
  levs <- levels(data[, variable])
  plotdata <- vector(mode="list", length=length(levs))
  for (i in seq_len(length(levs))) {
    cif_lev <- pseudo[data[, variable]==levs[i], ]
    cif_lev <- apply(cif_lev, 2, stats::weighted.mean,
                     w=weights[data[, variable]==levs[i]],
                     na.rm=TRUE)

    data_temp <- data.frame(time=times, cif=cif_lev, group=levs[i])

    if (conf_int) {
      # approximate variance by calculating a
      # weighted version of the variance
      cif_lev <- pseudo[data[, variable]==levs[i], ]

      cif_sd <- apply(cif_lev, 2, weighted.var.se,
                      w=weights[data[, variable]==levs[i]],
                      na.rm=TRUE, se_method=se_method)
      data_temp$se <- sqrt(cif_sd)

      cif_cis <- confint_surv(surv=data_temp$cif, se=data_temp$se,
                              conf_level=conf_level, conf_type="plain")
      data_temp$ci_lower <- cif_cis$left
      data_temp$ci_upper <- cif_cis$right

    }

    plotdata[[i]] <- data_temp
  }
  plotdata <- as.data.frame(dplyr::bind_rows(plotdata))
  rownames(plotdata) <- NULL

  output <- list(plotdata=plotdata,
                 pseudo_values=pseudo,
                 weights=weights)
  class(output) <- "adjustedcif.method"

  return(output)
}

## Direct Adjustment
#' @export
cif_direct <- function(data, variable, ev_time, event, cause, conf_int,
                       conf_level=0.95, times, outcome_model,
                       verbose=FALSE, predict_fun=NULL, ...) {

  # Using a Cause-Specific-Cox Model
  if (inherits(outcome_model, "CauseSpecificCox") & is.null(predict_fun)) {

    cif <- riskRegression::ate(event=outcome_model,
                               treatment=variable,
                               data=data,
                               estimator="Gformula",
                               times=times,
                               se=conf_int,
                               verbose=verbose,
                               cause=cause,
                               ...)
    plotdata <- data.frame(time=cif$meanRisk$time,
                           cif=cif$meanRisk$estimate,
                           group=cif$meanRisk$treatment)

    if (conf_int) {
      plotdata$se <- cif$meanRisk$se

      confint.ate <- utils::getFromNamespace("confint.ate", "riskRegression")

      cis <- confint.ate(cif, level=conf_level)$meanRisk
      plotdata$ci_lower <- cis$lower
      plotdata$ci_upper <- cis$upper
    }

    output <- list(plotdata=plotdata,
                   ate_object=cif)
    class(output) <- "adjustedcif.method"

  # Using a Fine & Gray Model
  } else {

    plotdata <- cif_g_comp(outcome_model=outcome_model,
                           data=data,
                           variable=variable,
                           times=times,
                           predict_fun=predict_fun,
                           cause=cause,
                           ...)

    output <- list(plotdata=plotdata)
    class(output) <- "adjustedcif.method"

  }

  return(output)
}

## using models other than coxph in method="direct" for
## survival endpoints
cif_g_comp <- function(outcome_model, data, variable, times,
                       predict_fun, cause, ...) {

  row_creation <- TRUE

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
                              cause=cause,
                              ...)
    # using predictRisk
    # NOTE: 'ranger' and 'fitSmoothHazard' accept a cause in predictRisk
    #       but don't actually work here.
    } else if (inherits(outcome_model, c("prodlim", "FGR", "rfsrc",
                                         "riskRegression", "ARR",
                                         "hal9001"))) {

      surv_lev <- riskRegression::predictRisk(object=outcome_model,
                                              newdata=data_temp,
                                              times=times,
                                              cause=cause,
                                              ...)

    # for fastCrr in fastcmprsk
    } else if (inherits(outcome_model, "fcrr")) {
      # get model matrix
      mod_vars <- all.vars(outcome_model$call[[2]])
      mod_form <- paste0(" ~ ", paste0(mod_vars, collapse=" + "))
      mod_data <- as.data.frame(stats::model.matrix(
        stats::as.formula(mod_form), data=data_temp))
      mod_data <- mod_data[, 2:ncol(mod_data)]

      # calculate average CIF
      rel_cols <- colnames(outcome_model$df[, 3:ncol(outcome_model$df)])
      surv_lev <- average_CIF_fccr(outcome_model, mod_data[, rel_cols])
      row <- data.frame(time=outcome_model$uftime,
                        cif=surv_lev,
                        group=levs[i])

      row_creation <- FALSE
    # for comp.risk in timereg
    } else if (inherits(outcome_model, "comprisk")) {
      surv_lev <- stats::predict(outcome_model, newdata=data_temp,
                                 times=times)$P1
    # using the S3 predict method
    } else {
      surv_lev <- tryCatch(
        expr={stats::predict(outcome_model,
                             newdata=data_temp,
                             times=times,
                             cause=cause,
                             ...)},
        error=function(e) {
          stop("The following error occured using",
               " the default S3 predict method: '", e,
               "' Specify a valid 'predict_fun' or",
               " use a different model. See details.")
          }
      )
    }

    # take arithmetic mean of predictions and add those to the
    # output object
    if (row_creation) {
      surv_lev <- apply(X=surv_lev, MARGIN=2, FUN=mean, na.rm=TRUE)
      row <- data.frame(time=times, cif=surv_lev, group=levs[i])
    }
    plotdata[[i]] <- row
  }
  plotdata <- as.data.frame(dplyr::bind_rows(plotdata))
  row.names(plotdata) <- seq_len(nrow(plotdata))

  return(plotdata)
}

## Matching
#' @export
cif_matching <- function(data, variable, ev_time, event, cause, conf_int,
                         conf_level=0.95, times, treatment_model,
                         gtol=0.001, ...) {

  # if it's a factor, turn it into numeric
  if (is.factor(data[, variable])) {
    levs <- levels(data[, variable])
    data[, variable] <- ifelse(data[, variable]==levs[1], 0, 1)
  } else {
    levs <- unique(data[, variable])
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
  m_dat <- rbind(data[rr$index.treated,], data[rr$index.control,])

  # estimate cif
  plotdata <- cif_aalen_johansen(data=m_dat,
                                 variable=variable,
                                 ev_time=ev_time,
                                 event=event,
                                 cause=cause,
                                 conf_int=FALSE,
                                 conf_level=conf_level)$plotdata

  if (!is.null(times)) {
    plotdata$se <- NULL
    plotdata <- specific_times(plotdata, times, cif=TRUE)
  }

  # get factor levels back
  plotdata$group[plotdata$group==0] <- levs[1]
  plotdata$group[plotdata$group==1] <- levs[2]

  output <- list(plotdata=plotdata,
                 match_object=rr)
  class(output) <- "adjustedcif.method"

  return(output)
}

## Using Augmented Inverse Probability of Treatment Weighting
#' @export
cif_aiptw <- function(data, variable, ev_time, event, cause, conf_int,
                      conf_level=0.95, times, outcome_model,
                      treatment_model, censoring_model=NULL,
                      verbose=FALSE, ...) {

  # defaults for input models
  if (is.null(censoring_model)) {
    form <- paste0("survival::Surv(", ev_time, ", ", event, "==0) ~ 1")
    censoring_model <- survival::coxph(stats::as.formula(form), data=data,
                                       x=TRUE, y=TRUE)
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
                               cause=cause)
  # calculate CIF estimates
  plotdata <- data.frame(time=curve$meanRisk$time,
                         cif=curve$meanRisk$estimate,
                         group=curve$meanRisk$treatment)
  if (conf_int) {
    plotdata$se <- curve$meanRisk$se

    confint.ate <- utils::getFromNamespace("confint.ate", "riskRegression")

    cis <- confint.ate(curve, level=conf_level, ci=TRUE)$meanRisk
    plotdata$ci_lower <- cis$lower
    plotdata$ci_upper <- cis$upper
  }

  output <- list(plotdata=plotdata,
                 ate_object=curve)
  class(output) <- "adjustedcif.method"

  return(output)
}


## Using Pseudo Observations and Direct Adjustment
#' @export
cif_direct_pseudo <- function(data, variable, ev_time, event, cause,
                              conf_int=FALSE, conf_level=0.95, times,
                              outcome_vars, type_time="factor", spline_df=5) {

  # estimate pseudo observations
  hist_formula <- stats::as.formula(paste("prodlim::Hist(", ev_time, ", ",
                                          event, ") ~ 1"))
  pseudo <- prodlim::jackknife(prodlim::prodlim(hist_formula, data=data),
                               times=times, cause=cause)

  # remove "variable" from outcome_vars because it is always included
  outcome_vars <- outcome_vars[outcome_vars!=variable]

  # some constants
  len <- length(times)
  n <- nrow(data)
  group <- data[, variable]

  # create data for geese
  Sdata <- data.frame(yi=1-c(pseudo),
                      group=rep(group, len),
                      vtime=rep(times, rep(n, len)),
                      id=rep(1:n, len))
  for (col in outcome_vars) {
    Sdata[,col] <- rep(data[,col], len)
  }

  if (type_time=="factor") {
    Sdata$vtime <- as.factor(Sdata$vtime)

    if (length(times)==1) {
      geese_formula <- paste("yi ~ ", paste(outcome_vars, collapse=" + "),
                             " + group")
    } else {
      geese_formula <- paste("yi ~ vtime + ", paste(outcome_vars,
                                                    collapse=" + "),
                             " + group")
  }

  } else if (type_time=="bs") {
    geese_formula <- paste("yi ~ splines::bs(vtime, df=", spline_df, ") + ",
                           paste(outcome_vars, collapse=" + "), " + group")
  } else if (type_time=="ns") {
    geese_formula <- paste("yi ~ splines::ns(vtime, df=", spline_df, ") + ",
                           paste(outcome_vars, collapse=" + "), " + group")
  }

  # call geese
  geese_mod <- geepack::geese(stats::as.formula(geese_formula), scale.fix=TRUE,
                              data=Sdata, family=gaussian, id=id, jack=FALSE,
                              mean.link="cloglog", corstr="independence")

  # initialize outcome df list
  levs <- levels(data[,variable])
  plotdata <- vector(mode="list", length=length(levs))

  # do direct adjustment
  for (i in seq_len(length(levs))) {

    Sdata$group <- factor(levs[i], levels=levs)
    pred <- geese_predictions(geese_mod, Sdata, times=times, n=n)

    m <- exp(-exp(pred))
    cif <- apply(m, 2, mean, na.rm=TRUE)

    plotdata[[i]] <- data.frame(time=times, cif=cif, group=levs[i])

  }
  plotdata <- dplyr::bind_rows(plotdata)
  rownames(plotdata) <- NULL

  output <- list(plotdata=plotdata,
                 pseudo_values=pseudo,
                 geese_model=geese_mod)
  class(output) <- "adjustedcif.method"

  return(output)
}

## Using AIPTW with Pseudo Observations
#' @export
cif_aiptw_pseudo <- function(data, variable, ev_time, event, cause,
                             conf_int, conf_level=0.95, times,
                             outcome_vars, treatment_model,
                             type_time="factor", spline_df=5) {
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
  hist_formula <- stats::as.formula(paste("prodlim::Hist(", ev_time, ", ",
                                          event, ") ~ 1"))
  pseudo <- prodlim::jackknife(prodlim::prodlim(hist_formula, data=data),
                               times=times, cause=cause)
  # create data for geese
  Sdata <- data.frame(yi=c(pseudo),
                      group=rep(group, len),
                      vtime=rep(times, rep(n, len)),
                      id=rep(1:n, len))

  outcome_vars <- outcome_vars[outcome_vars != variable]
  for (col in outcome_vars) {
    Sdata[,col] <- rep(data[, col], len)
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

  # call geese
  geese_mod <- geepack::geese(stats::as.formula(geese_formula), scale.fix=TRUE,
                              data=Sdata, family=gaussian, id=id, jack=FALSE,
                              mean.link="cloglog", corstr="independence")

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

    cif <- apply(dr, 2, mean, na.rm=TRUE)

    if (conf_int) {

      pseudo_dr_se <- function(x, n, na.rm) {
        sqrt(stats::var(x, na.rm=na.rm) / n)
      }
      cif_se <- apply(dr, 2, pseudo_dr_se, n=n, na.rm=TRUE)

      cif_ci <- confint_surv(surv=cif, se=cif_se, conf_level=conf_level,
                             conf_type="plain")

      plotdata[[i]] <- data.frame(time=times, cif=cif, group=levs[i],
                                  se=cif_se, ci_lower=cif_ci$left,
                                  ci_upper=cif_ci$right)

    } else {
      plotdata[[i]] <- data.frame(time=times, cif=cif, group=levs[i])
    }

  }
  plotdata <- dplyr::bind_rows(plotdata)
  rownames(plotdata) <- NULL

  output <- list(plotdata=plotdata,
                 pseudo_values=pseudo,
                 geese_model=geese_mod)
  class(output) <- "adjustedcif.method"

  return(output)
}

## Targeted Maximum Likelihood Estimation
#' @export
cif_tmle <- function(data, variable, ev_time, event, cause, conf_int,
                     conf_level=0.95, times, adjust_vars=NULL,
                     SL.ftime=NULL, SL.ctime=NULL, SL.trt=NULL,
                     glm.ftime=NULL, glm.ctime=NULL, glm.trt=NULL,
                            ...) {
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

  cif_0 <- unlist(lapply(tpfit, function(x) {x$est[1]}))
  cif_1 <- unlist(lapply(tpfit, function(x) {x$est[2]}))

  # put together
  plotdata <- data.frame(time=rep(times, 2),
                         cif=c(cif_0, cif_1),
                         group=c(rep(levs[1], length(times)),
                                 rep(levs[2], length(times))))

  if (conf_int) {

    var_0 <- unlist(lapply(tpfit, function(x) {x$var[1]}))
    var_1 <- unlist(lapply(tpfit, function(x) {x$var[2]}))

    plotdata$se <- sqrt(c(var_0, var_1))

    confint.tp.survtmle <- utils::getFromNamespace("confint.tp.survtmle",
                                                   "survtmle")
    survtmle_ci <- confint.tp.survtmle(tpfit, level=conf_level)

    plotdata$ci_lower <- c(survtmle_ci$`0 1`[,1], survtmle_ci$`1 1`[,1])
    plotdata$ci_upper <- c(survtmle_ci$`0 1`[,2], survtmle_ci$`1 1`[,2])

  }

  output <- list(plotdata=plotdata,
                 survtmle_object=fit,
                 survtmle.timepoints_object=tpfit)
  class(output) <- "adjustedcif.method"

  return(output)
}
