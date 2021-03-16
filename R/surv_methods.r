# Assign to global, to get rid off devtools::check() note
utils::globalVariables(c("gaussian", "id"))

## simple Kaplan-Meier estimate
#' @export
surv_method_km <- function(data, variable, ev_time, event, conf_int,
                           conf_level=0.95, ...) {

  form <- paste0("survival::Surv(", ev_time, ", ", event, ") ~ ", variable)

  surv <- survival::survfit(stats::as.formula(form), data=data, se.fit=conf_int,
                            conf.int=conf_level, ...)
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

  return(plotdata)
}

## IPTW Kaplan-Meier estimate
# TODO: - standard deviation calculation is a little off,
#       - CI seem fine for group==1 but way too small otherwise
#' @export
surv_method_iptw_km <- function(data, variable, ev_time, event, conf_int,
                                conf_level=0.95, treatment_model,
                                weight_method="ps", stabilize=T, ...) {

  # get weights
  if (is.numeric(treatment_model)) {
    weights <- treatment_model
  } else {
    weights <- get_iptw_weights(data=data, treatment_model=treatment_model,
                                weight_method=weight_method,
                                variable=variable, stabilize=stabilize, ...)
  }

  # weighted Kaplan-Meier estimator
  levs <- levels(data[,variable])
  plotdata <- vector(mode="list", length=length(levs))

  for (i in 1:length(levs)) {

    data_lev <- data[data[,variable]==levs[i],]
    weights_lev <- weights[data[,variable]==levs[i]]

    times <- sort(unique(data_lev[,ev_time][data_lev[,event]==1]))
    max_t <- max(data_lev[,ev_time])

    if (!0 %in% times) {
      times <- c(0, times)
    }
    if (!max_t %in% times) {
      times <- c(times, max_t)
    }

    d_j <- sapply(times, function(x){sum(weights_lev[data_lev[,ev_time]==x &
                                                   data_lev[,event]==1])})
    Y_j <- sapply(times, function(x){sum(weights_lev[data_lev[,ev_time]>=x])})
    S_t <- cumprod((Y_j-d_j)/Y_j)

    plotdata[[i]] <- data.frame(time=times, surv=S_t, group=levs[i],
                                d_j=d_j, Y_j=Y_j)
  }
  plotdata <- dplyr::bind_rows(plotdata)

  # variance, confidence interval
  if (conf_int) {
    plotdata$sd <- NA
    plotdata$s_j <- 1 - (plotdata$d_j / plotdata$Y_j)

    # get propensity scores
    if (inherits(treatment_model, "glm")) {
      ps_score <- treatment_model$fitted.values
    } else if (inherits(treatment_model, "multinom")) {
      predict.multinom <- utils::getFromNamespace("predict.multinom", "nnet")
      ps_score <- predict.multinom(treatment_model, newdata=data, type="probs")
    }

    for (i in 1:length(levs)) {

      # relevant propensity scores
      if (length(levs) > 2) {
        ps_score_lev <- ps_score[,levs[i]]
      } else if (i==1) {
        ps_score_lev <- 1 - ps_score
      } else if (i==2) {
        ps_score_lev <- ps_score
      }

      adj_km_lev <- plotdata[plotdata$group==levs[i],]

      # calculate Mj for every relevant point in time
      adj_km_lev$Mj <- sapply(adj_km_lev$time, calc_Mj, data=data,
                              ev_time=ev_time, ps_score=ps_score_lev)

      # calculate variance at each point in time
      iptw_km_var <- sapply(adj_km_lev$time, calc_iptw_km_var,
                            adj_km=adj_km_lev)
      plotdata$sd[plotdata$group==levs[i]] <- sqrt(iptw_km_var)
    }
    plotdata$s_j <- NULL

    surv_cis <- confint_surv(surv=plotdata$surv, se=plotdata$sd,
                             conf_level=conf_level, conf_type="plain")
    plotdata$ci_lower <- surv_cis$left
    plotdata$ci_upper <- surv_cis$right
  }

  plotdata$d_j <- NULL
  plotdata$Y_j <- NULL

  return(plotdata)
}

## IPTW with univariate cox-model
#' @export
surv_method_iptw_cox <- function(data, variable, ev_time, event, conf_int,
                                 conf_level=0.95, treatment_model,
                                 weight_method="ps", stabilize=T, ...) {

  # get weights
  if (is.numeric(treatment_model)) {
    weights <- treatment_model
  } else {
    weights <- get_iptw_weights(data=data, treatment_model=treatment_model,
                                weight_method=weight_method,
                                variable=variable, stabilize=stabilize, ...)
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

  return(plotdata)
}

## Using Pseudo-Observations and IPTW
#' @export
surv_method_iptw_pseudo <- function(data, variable, ev_time, event, conf_int,
                                    conf_level=0.95, times, treatment_model,
                                    weight_method="ps", stabilize=T,
                                    se_method="cochrane", ...) {
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
                               times=times)

  # take weighted mean
  levs <- levels(data[,variable])
  plotdata <- vector(mode="list", length=length(levs))
  for (i in 1:length(levs)) {
    surv_lev <- pseudo[data[,variable]==levs[i],]
    surv_lev <- apply(surv_lev, 2, stats::weighted.mean,
                      w=weights[data[,variable]==levs[i]],
                      na.rm=T)

    data_temp <- data.frame(time=times, surv=surv_lev, group=levs[i])

    if (conf_int) {
      # approximate variance by calculating a
      # weighted version of the variance
      surv_lev <- pseudo[data[,variable]==levs[i],]

      surv_sd <- apply(surv_lev, 2, weighted.var.se,
                       w=weights[data[,variable]==levs[i]],
                       na.rm=T, se_method=se_method)
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

  return(plotdata)
}

## Direct Adjustment
#' @export
surv_method_direct <- function(data, variable, ev_time, event, conf_int,
                               conf_level=0.95, times, outcome_model,
                               verbose=F, ...) {

  surv <- riskRegression::ate(event=outcome_model, treatment=variable,
                              data=data, estimator="Gformula",
                              times=times, se=conf_int, verbose=verbose, ...)
  plotdata <- data.frame(time=surv$meanRisk$time,
                         surv=1 - surv$meanRisk$estimate,
                         group=surv$meanRisk$treatment)

  if (conf_int) {
    plotdata$se <- surv$meanRisk$se

    confint.ate <- utils::getFromNamespace("confint.ate", "riskRegression")

    cis <- confint.ate(surv, level=conf_level)$meanRisk
    plotdata$ci_lower <- 1 - cis$lower
    plotdata$ci_upper <- 1 - cis$upper
  }

  return(plotdata)
}

## Using propensity score matching
# TODO: variance calculation is off
#' @export
surv_method_matching <- function(data, variable, ev_time, event, conf_int,
                                 conf_level=0.95, treatment_model,
                                 stabilize=T, ...) {

  if (is.numeric(treatment_model)) {
    ps_score <- treatment_model
  } else {
    ps_score <- stats::predict.glm(treatment_model, newdata=data,
                                   type="response")
  }

  rr <- Matching::Match(Tr=data[, variable], X=ps_score, estimand="ATE", ...)
  m_dat <- rbind(data[rr$index.treated,], data[rr$index.control,])

  weights <- rr$weights

  if (stabilize) {
    weights <- weights * length(weights) / sum(weights)
  }
  m_dat$match_weights <- c(weights, weights)

  # estimate survival curve
  form <- paste0("survival::Surv(", ev_time, ", ", event, ") ~ ", variable)
  surv <- survival::survfit(stats::as.formula(form), data=m_dat, se.fit=conf_int,
                            conf.int=conf_level, weights=m_dat$match_weights,
                            robust=T)
  plotdata <- data.frame(time=surv$time,
                         surv=surv$surv,
                         group=c(rep(0, surv$strata[1]),
                                 rep(1, surv$strata[2])))

  if (conf_int) {
    plotdata$se <- surv$std.err
    plotdata$ci_lower <- surv$lower
    plotdata$ci_upper <- surv$upper
  }

  return(plotdata)
}

## Using Augmented Inverse Probability of Treatment Weighting
#' @export
surv_method_aiptw <- function(data, variable, ev_time, event, conf_int,
                              conf_level=0.95, times, outcome_model=NULL,
                              treatment_model=NULL, censoring_model=NULL,
                              verbose=F, ...) {

  # defaults for input models
  if (is.null(censoring_model)) {
    form <- paste0("survival::Surv(", ev_time, ", ", event, "==0) ~ 1")
    censoring_model <- survival::coxph(stats::as.formula(form), data=data, x=T)
  }
  if (is.null(treatment_model)) {
    form <- paste0(variable, " ~ 1")
    treatment_model <- stats::glm(stats::as.formula(form), data=data,
                                  family="binomial")
  }
  if (is.null(outcome_model)) {
    form <- paste0("survival::Surv(", ev_time, ", ", event, ") ~ 1")
    outcome_model <- survival::coxph(stats::as.formula(form), data=data, x=T)
  }

  # estimate AIPTW cumulative incidence
  curve <- riskRegression::ate(event=outcome_model,
                               treatment=treatment_model,
                               censor=censoring_model,
                               data=data,
                               times=times,
                               se=conf_int,
                               verbose=verbose,
                               estimator="AIPTW,AIPCW")
  # calculate survival estimates
  plotdata <- data.frame(time=curve$meanRisk$time,
                         surv=1 - curve$meanRisk$estimate,
                         group=curve$meanRisk$treatment)
  if (conf_int) {
    plotdata$se <- curve$meanRisk$se

    confint.ate <- utils::getFromNamespace("confint.ate", "riskRegression")

    cis <- confint.ate(curve, level=conf_level, ci=T)$meanRisk
    plotdata$ci_lower <- 1 - cis$lower
    plotdata$ci_upper <- 1 - cis$upper
  }

  return(plotdata)
}

## Using Pseudo Observations and Direct Adjustment
# TODO: fails when len==1: either add linear regression or stop()
#' @export
surv_method_direct_pseudo <- function(data, variable, ev_time, event, times,
                                      outcome_vars, type_time="factor",
                                      spline_df=10) {
  # some constants
  len <- length(times)
  n <- nrow(data)
  group <- data[,variable]

  # estimate pseudo observations
  hist_formula <- stats::as.formula(paste("prodlim::Hist(", ev_time, ", ",
                                    event, ") ~ 1"))
  pseudo <- prodlim::jackknife(prodlim::prodlim(hist_formula, data=data),
                               times=times)
  # create data for geese
  Sdata <- data.frame(yi=c(pseudo),
                      group=rep(group, len),
                      vtime=rep(times, rep(n, len)),
                      id=rep(1:n, len))
  for (col in outcome_vars) {
    Sdata[,col] <- rep(data[,col], len)
  }

  Sdata <- Sdata[,!is.na(Sdata$yi)]

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
    surv <- apply(m, 2, mean, na.rm=T)

    plotdata[[i]] <- data.frame(time=times, surv=surv, group=levs[i])

  }
  plotdata <- dplyr::bind_rows(plotdata)
  rownames(plotdata) <- NULL

  return(plotdata)
}

## Using Pseudo Observations Doubly-Robust, Wang (2018)
#' @export
surv_method_aiptw_pseudo <- function(data, variable, ev_time, event, conf_int,
                                     conf_level=0.95, times, outcome_vars,
                                     treatment_model, type_time="factor",
                                     spline_df=10) {
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
                               times=times)
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

    surv <- apply(dr, 2, mean, na.rm=T)

    if (conf_int) {

      pseudo_dr_se <- function(x, n, na.rm) {
        sqrt(stats::var(x, na.rm=na.rm) / n)
      }
      surv_se <- apply(dr, 2, pseudo_dr_se, n=n, na.rm=T)

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

  return(plotdata)
}

## Using Empirical Likelihood Estimation
#' @export
surv_method_el <- function(data, variable, ev_time, event,
                           times, treatment_vars, moment="first",
                           standardize=F, ...) {

  el_0 <- adjKMtest::el.est(y=data[, ev_time],
                            delta=data[, event],
                            treat=data[, variable],
                            x=as.matrix(data[, treatment_vars]),
                            treat.select=0,
                            t=times,
                            psix_moment=moment,
                            standardize=standardize,
                            get.sd=F,
                            ...)
  el_1 <- adjKMtest::el.est(y=data[, ev_time],
                            delta=data[, event],
                            treat=data[, variable],
                            x=as.matrix(data[, treatment_vars]),
                            treat.select=1,
                            t=times,
                            psix_moment=moment,
                            standardize=standardize,
                            get.sd=F,
                            ...)

  plotdata <- data.frame(time=c(times, times),
                         surv=c(el_0$St, el_1$St),
                         group=c(rep(0, length(el_0$St)),
                                 rep(1, length(el_0$St))))

  return(plotdata)
}

## Targeted Maximum Likelihood Estimation
# TODO: Fails with multiple time points cause of weird evaluation in "timepoints"
#' @export
surv_method_tmle <- function(data, variable, ev_time, event, conf_int,
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

  s_0 <- 1 - unlist(lapply(tpfit, function(x) {x$est[1]}))
  s_1 <- 1 - unlist(lapply(tpfit, function(x) {x$est[2]}))

  # put together
  plotdata <- data.frame(time=rep(times, 2),
                         surv=c(s_0, s_1),
                         group=c(rep(0, length(times)),
                                 rep(1, length(times))))

  if (conf_int) {

    var_0 <- unlist(lapply(tpfit, function(x) {x$var[1]}))
    var_1 <- unlist(lapply(tpfit, function(x) {x$var[2]}))

    plotdata$sd <- sqrt(c(var_0, var_1))

    confint.tp.survtmle <- utils::getFromNamespace("confint.tp.survtmle",
                                                   "survtmle")
    survtmle_ci <- confint.tp.survtmle(tpfit, level=conf_level)

    plotdata$ci_lower <- 1 - c(survtmle_ci$`0 1`[,1], survtmle_ci$`1 1`[,1])
    plotdata$ci_upper <- 1 - c(survtmle_ci$`0 1`[,2], survtmle_ci$`1 1`[,2])

  }

  return(plotdata)
}

## One-Step Targeted Maximum Likelihood Estimation
#' @export
surv_method_ostmle <- function(data, variable, ev_time, event, conf_int,
                               conf_level=0.95, times, adjust_vars=NULL,
                               SL.ftime=NULL, SL.ctime=NULL, SL.trt=NULL,
                               epsilon=1, max_num_iteration=100,
                               psi_moss_method="l2", tmle_tolerance=NULL,
                               gtol=1e-3) {
  if (is.null(adjust_vars)) {
    all_covars <- colnames(data)
    all_covars <- all_covars[!all_covars %in% c(variable, ev_time, event)]
    adjust_vars <- data[, all_covars]
  } else {
    adjust_vars <- data[, adjust_vars]
  }

  # time point grid
  k_grid <- 1:max(data[,ev_time])

  # create initial fit object
  sl_fit <- MOSS::initial_sl_fit(
    T_tilde=data[, ev_time],
    Delta=data[, event],
    A=data[, variable],
    W=adjust_vars,
    t_max=max(data[,ev_time]),
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
  moss_hazard_fit_1 <- MOSS::MOSS_hazard$new(
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
  moss_hazard_fit_0 <- MOSS::MOSS_hazard$new(
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
                         group=c(rep(0, length(k_grid)),
                                 rep(1, length(k_grid))))
  # keep only time points in times
  plotdata <- plotdata[which(plotdata$time %in% times),]

  if (conf_int) {

    # for treatment == 1
    s_1 <- as.vector(psi_moss_hazard_1)
    eic_fit <- MOSS::eic$new(
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
    eic_fit <- MOSS::eic$new(
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

  return(plotdata)
}
