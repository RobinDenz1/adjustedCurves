##################### Similar Methods as in Survival Context ###################

# A tradeoff had to be made here. As can clearly be seen, some cif_method
# functions are nearly identical to the surv_method equivalents. Instead
# of writing one more general function, we choose to write separate functions
# to make the documentation and usage easier, while making the code itself
# a little messier with unnecessary repetition.

## Aalen-Johansen estimator
#' @export
cif_method_aalen <- function(data, variable, ev_time, event, cause, conf_int,
                             conf_level=0.95, ...) {

  cif <- cmprsk::cuminc(ftime=data[,ev_time],
                        fstatus=data[,event],
                        group=data[,variable],
                        ...)

  levs <- unique(data[,variable])
  cif_names <- paste(levs, cause)
  plotdata <- vector(mode="list", length=length(cif_names))
  for (i in 1:length(cif_names)) {
    plotdata[[i]] <- data.frame(time=cif[[cif_names[i]]]$time,
                                cif=cif[[cif_names[i]]]$est,
                                group=levs[i],
                                se=cif[[cif_names[i]]]$var)
  }
  plotdata <- as.data.frame(dplyr::bind_rows(plotdata))

  if (conf_int) {
    cif_cis <- confint_surv(surv=plotdata$cif, se=plotdata$se,
                             conf_level=conf_level, conf_type="plain")
    plotdata$ci_lower <- cif_cis$left
    plotdata$ci_upper <- cif_cis$right
  }

  return(plotdata)
}

## IPTW
#' @export
cif_method_iptw <- function(data, variable, ev_time, event, cause, conf_int,
                            conf_level=0.95, times, treatment_model,
                            censoring_model=NULL, verbose=F, ...) {
  # empty censoring model if not specified
  if (is.null(censoring_model)) {
    form <- paste0("survival::Surv(", ev_time, ", ", event, "==0) ~ 1")
    censoring_model <- survival::coxph(stats::as.formula(form), data=data,
                                       x=T, y=T)
  }

  cif <- riskRegression::ate(event=c(ev_time, event), treatment=treatment_model,
                             data=data, estimator="IPTW",
                             times=times, se=conf_int, verbose=verbose,
                             cause=cause, censor=censoring_model, ...)
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

  return(plotdata)
}

# IPTW pseudo
#' @export
cif_method_iptw_pseudo <- function(data, variable, ev_time, event, cause,
                                   conf_int, conf_level=0.95, times,
                                   treatment_model, weight_method="ps",
                                   stabilize=T, se_method="cochrane", ...) {
  # get weights
  if (is.numeric(treatment_model)) {
    weights <- treatment_model
  } else {
    weights <- get_iptw_weights(data=data, treatment_model=treatment_model,
                                weight_method=weight_method,
                                variable=variable, stabilize=stabilize, ...)
  }

  # estimate pseudo observations
  hist_formula <- stats::as.formula(paste("prodlim::Hist(", ev_time, ", ",
                                          event, ") ~ 1"))
  pseudo <- prodlim::jackknife(prodlim::prodlim(hist_formula, data=data),
                               times=times, cause=cause)

  # take weighted mean
  levs <- levels(data[,variable])
  plotdata <- vector(mode="list", length=length(levs))
  for (i in 1:length(levs)) {
    cif_lev <- pseudo[data[,variable]==levs[i],]
    cif_lev <- apply(cif_lev, 2, stats::weighted.mean,
                     w=weights[data[,variable]==levs[i]],
                     na.rm=T)

    data_temp <- data.frame(time=times, cif=cif_lev, group=levs[i])

    if (conf_int) {
      # approximate variance by calculating a
      # weighted version of the variance
      cif_lev <- pseudo[data[,variable]==levs[i],]

      cif_sd <- apply(cif_lev, 2, weighted.var.se,
                      w=weights[data[,variable]==levs[i]],
                      na.rm=T, se_method=se_method)
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

  return(plotdata)
}

## Direct Adjustment
#' @export
cif_method_direct <- function(data, variable, ev_time, event, cause, conf_int,
                              conf_level=0.95, times, outcome_model,
                              verbose=F, ...) {

  cif <- riskRegression::ate(event=outcome_model, treatment=variable,
                             data=data, estimator="Gformula",
                             times=times, se=conf_int, verbose=verbose,
                             cause=cause, ...)
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

  return(plotdata)
}

## Matching
# TODO: variance calculation is off
#' @export
cif_method_matching <- function(data, variable, ev_time, event, cause, conf_int,
                                conf_level=0.95, treatment_model, ...) {

  if (is.numeric(treatment_model)) {
    ps_score <- treatment_model
  } else {
    ps_score <- stats::predict.glm(treatment_model, newdata=data,
                                   type="response")
  }

  rr <- Matching::Match(Tr=data[, variable], X=ps_score, estimand="ATE", ...)
  m_dat <- rbind(data[rr$index.treated,], data[rr$index.control,])

  # estimate cif
  plotdata <- cif_method_aalen(data=m_dat, variable=variable,
                               ev_time=ev_time, event=event,
                               cause=cause, conf_int=conf_int,
                               conf_level=conf_level)

  return(plotdata)
}

## Using Augmented Inverse Probability of Treatment Weighting
#' @export
cif_method_aiptw <- function(data, variable, ev_time, event, cause, conf_int,
                             conf_level=0.95, times, outcome_model=NULL,
                             treatment_model=NULL, censoring_model=NULL,
                             verbose=F, ...) {

  # defaults for input models
  if (is.null(censoring_model)) {
    form <- paste0("survival::Surv(", ev_time, ", ", event, "==0) ~ 1")
    censoring_model <- survival::coxph(stats::as.formula(form), data=data,
                                       x=T, y=T)
  }
  if (is.null(treatment_model)) {
    form <- paste0(variable, " ~ 1")
    treatment_model <- stats::glm(stats::as.formula(form), data=data,
                                  family="binomial")
  }
  if (is.null(outcome_model)) {
    form <- paste0("prodlim::Hist(", ev_time, ", ", event, ") ~ 1")
    outcome_model <- riskRegression::CSC(stats::as.formula(form), data=data, x=T)
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

    cis <- confint.ate(curve, level=conf_level, ci=T)$meanRisk
    plotdata$ci_lower <- cis$lower
    plotdata$ci_upper <- cis$upper
  }

  return(plotdata)
}


## Using Pseudo Observations and Direct Adjustment
# TODO: fails when len==1: either add linear regression or stop()
#' @export
cif_method_direct_pseudo <- function(data, variable, ev_time, event, cause,
                                     times, outcome_vars, type_time="factor",
                                     spline_df=10) {
  # some constants
  len <- length(times)
  n <- nrow(data)
  group <- data[,variable]

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
  for (col in outcome_vars) {
    Sdata[,col] <- rep(data[,col], len)
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
                              data=Sdata, family=gaussian, id=id, jack=F,
                              mean.link="cloglog", corstr="independence")

  # initialize outcome df list
  levs <- levels(data[,variable])
  plotdata <- vector(mode="list", length=length(levs))

  # do direct adjustment
  for (i in 1:length(levs)) {

    Sdata$group <- factor(levs[i], levels=levs)
    pred <- geese_predictions(geese_mod, Sdata, times=times, n=n)

    m <- 1 - exp(-exp(pred))
    cif <- apply(m, 2, mean, na.rm=T)

    plotdata[[i]] <- data.frame(time=times, cif=cif, group=levs[i])

  }
  plotdata <- dplyr::bind_rows(plotdata)
  rownames(plotdata) <- NULL

  return(plotdata)
}

## Using AIPTW with Pseudo Observations
#' @export
cif_method_aiptw_pseudo <- function(data, variable, ev_time, event, cause,
                                    conf_int, conf_level=0.95, times,
                                    outcome_vars, treatment_model,
                                    type_time="factor", spline_df=10) {
  # some constants
  len <- length(times)
  n <- nrow(data)
  group <- data[,variable]

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
    Sdata[,col] <- rep(data[,col], len)
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
                              data=Sdata, family=gaussian, id=id, jack=F,
                              mean.link="cloglog", corstr="independence")

  # get direct adjustment estimates
  levs <- levels(data[,variable])
  plotdata <- vector(mode="list", length=length(levs))

  for (i in 1:length(levs)) {

    Sdata$group <- factor(levs[i], levels=levs)
    pred <- geese_predictions(geese_mod, Sdata, times=times, n=n)

    m <- 1 - exp(-exp(pred))

    # augment estimates using propensity score
    if (length(levs) > 2) {
      group_ind <- ifelse(group==levs[i], 1, 0)
      ps_score_lev <- ps_score[,levs[i]]

      dr <- (pseudo*group_ind-(group_ind-ps_score_lev)*m)/ps_score_lev
      # if binary, use equation from the paper directly
    } else if (i == 1) {
      group <- ifelse(data[,variable]==levs[1], 0, 1)
      dr <- (pseudo*(1-group)+(group-ps_score)*m)/(1-ps_score)
    } else if (i == 2) {
      group <- ifelse(data[,variable]==levs[1], 0, 1)
      dr <- (pseudo*group-(group-ps_score)*m)/ps_score
    }

    cif <- apply(dr, 2, mean, na.rm=T)

    if (conf_int) {

      pseudo_dr_se <- function(x, n, na.rm) {
        sqrt(stats::var(x, na.rm=na.rm) / n)
      }
      cif_se <- apply(dr, 2, pseudo_dr_se, n=n, na.rm=T)

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

  return(plotdata)
}

## Targeted Maximum Likelihood Estimation
# TODO: Fails with multiple time points cause of weird evaluation in "timepoints"
#' @export
cif_method_tmle <- function(data, variable, ev_time, event, cause, conf_int,
                            conf_level=0.95, times, adjust_vars=NULL,
                            SL.ftime=NULL, SL.ctime=NULL, SL.trt=NULL,
                            glm.ftime=NULL, glm.ctime=NULL, glm.trt=NULL,
                            ...) {
  if (is.null(adjust_vars)) {
    all_covars <- colnames(data)
    all_covars <- all_covars[!all_covars %in% c(variable, ev_time, event)]
    adjust_vars <- data[,all_covars]
  } else {
    adjust_vars <- data[,adjust_vars]
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
    survtmle::timepoints(fit, times=times)
  }, warning=function(w) {
    if (startsWith(conditionMessage(w), "Using formula(x) is deprecated"))
      invokeRestart("muffleWarning")
  })

  cif_0 <- unlist(lapply(tpfit, function(x) {x$est[1]}))
  cif_1 <- unlist(lapply(tpfit, function(x) {x$est[2]}))

  # put together
  plotdata <- data.frame(time=rep(times, 2),
                         cif=c(cif_0, cif_1),
                         group=c(rep(0, length(times)),
                                 rep(1, length(times))))

  if (conf_int) {

    var_0 <- unlist(lapply(tpfit, function(x) {x$var[1]}))
    var_1 <- unlist(lapply(tpfit, function(x) {x$var[2]}))

    plotdata$sd <- sqrt(c(var_0, var_1))

    confint.tp.survtmle <- utils::getFromNamespace("confint.tp.survtmle",
                                                   "survtmle")
    survtmle_ci <- confint.tp.survtmle(tpfit, level=conf_level)

    plotdata$ci_lower <- c(survtmle_ci$`0 1`[,1], survtmle_ci$`1 1`[,1])
    plotdata$ci_upper <- c(survtmle_ci$`0 1`[,2], survtmle_ci$`1 1`[,2])

  }

  return(plotdata)
}
