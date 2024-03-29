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

## Direct Adjustment
#' @export
surv_direct <- function(data, variable, ev_time, event, conf_int,
                        conf_level=0.95, times, outcome_model,
                        verbose=FALSE, predict_fun=NULL, ...) {

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

      surv_lev <- quiet(riskRegression::predictRisk(object=outcome_model,
                                              newdata=data_temp,
                                              times=times,
                                              cause=1,
                                              ...))
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

  # Using a Fine & Gray or other Model
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

## using models other than CauseSpecificCox in method="direct" for
## competing-risks endpoints
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

      surv_lev <- quiet(riskRegression::predictRisk(object=outcome_model,
                                              newdata=data_temp,
                                              times=times,
                                              cause=cause,
                                              ...))

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

#### For Adjusted CIFs ####

## fastcmprsk

# multiply betas with covariates and take sum
apply_betas <- function(x, betas) {
  return(sum(x * betas))
}

# using the same method as in predict.fccr, but vectorised to get the
# predictions in matrix form first
average_CIF_fccr <- function(object, newdata) {

  # apply betas
  eff <- apply(X=newdata, MARGIN=1, FUN=apply_betas, betas=object$coef)

  # same as in predict.fccr
  CIF.hat <- t(outer(X=exp(eff), Y=object$breslowJump[, 2], FUN="*"))
  CIF.hat <- apply(X=CIF.hat, MARGIN=2, FUN=cumsum)
  CIF.hat <- apply(X=CIF.hat, MARGIN=2, FUN=function(x){1 - exp(-x)})

  # take average at each t
  CIF_z <- apply(X=CIF.hat, MARGIN=1, FUN=mean)

  return(CIF_z)
}

##
