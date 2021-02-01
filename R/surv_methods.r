## simple Kaplan-Meier estimate
#' @export
surv_method_km <- function(data, variable, ev_time, event, sd, ...) {

  form <- paste0("Surv(", ev_time, ", ", event, ") ~ ", variable)

  surv <- survfit(as.formula(form), data=data, se.fit=sd, ...)
  plotdata <- data.frame(time=surv$time,
                         surv=surv$surv,
                         group=c(rep(0, surv$strata[1]),
                                 rep(1, surv$strata[2])))
  if (sd) {
    plotdata$sd <- surv$std.err * sqrt(nrow(data))
  }

  return(plotdata)
}


## IPTW Kaplan-Meier estimate
#' @export
surv_method_iptw_km <- function(data, variable, ev_time, event, sd,
                                treatment_model, weight_method="ps", ...) {

  # get weights
  if (inherits(treatment_model, "formula") | inherits(treatment_model, "glm")) {
    weights <- get_iptw_weights(data, treatment_model, weight_method, ...)
  } else if (is.numeric(treatment_model)) {
    weights <- treatment_model
  }

  # weighted Kaplan-Meier estimator
  adj_km <- RISCA::ipw.survival(data[, ev_time], data[, event],
                                data[, variable], weights=weights)
  plotdata <- data.frame(time=adj_km[[1]]$times,
                         surv=adj_km[[1]]$survival,
                         group=adj_km[[1]]$variable)
  return(plotdata)
}

## IPTW with univariate cox-model
#' @export
surv_method_iptw_cox <- function(data, variable, ev_time, event, sd,
                                 treatment_model, weight_method="ps", ...) {

  # get weights
  if (inherits(treatment_model, "formula") | inherits(treatment_model, "glm")) {
    weights <- get_iptw_weights(data, treatment_model, weight_method, ...)
  } else if (is.numeric(treatment_model)) {
    weights <- treatment_model
  }

  # univariate, weighted cox model
  form <- paste0("Surv(", ev_time, ", ", event, ") ~ strata(", variable, ")")
  model <- survival::coxph(as.formula(form), weights=weights, data=data)
  surv <- survival::survfit(model, se.fit=sd)

  plotdata <- data.frame(time=surv$time,
                         surv=surv$surv,
                         group=c(rep(0, surv$strata[1]),
                                 rep(1, surv$strata[2])))

  if (sd) {
    plotdata$sd <- surv$std.err * sqrt(nrow(data))
  }

  return(plotdata)
}

## Using Pseudo-Observations and IPTW
#' @export
surv_method_iptw_pseudo <- function(data, variable, ev_time, event, sd,
                                    times, treatment_model,
                                    weight_method="ps", na.rm=F, ...) {
  # get weights
  if (inherits(treatment_model, "formula") | inherits(treatment_model, "glm")) {
    weights <- get_iptw_weights(data, treatment_model, weight_method, ...)
  } else if (is.numeric(treatment_model)) {
    weights <- treatment_model
  }

  # estimate pseudo observations
  hist_formula <- as.formula(paste("Hist(", ev_time, ", ", event, ") ~ 1"))
  pseudo <- prodlim::jackknife(prodlim(hist_formula, data=data),
                               times=times$overall)

  # take weighted mean
  surv_A <- pseudo[data[,variable]==0,]
  surv_B <- pseudo[data[,variable]==1,]
  surv_A <- apply(surv_A, 2, weighted.mean, w=weights[data[,variable]==0],
                  na.rm=na.rm)
  surv_B <- apply(surv_B, 2, weighted.mean, w=weights[data[,variable]==1],
                  na.rm=na.rm)
  # put together
  plotdata <- data.frame(time=rep(times$overall, 2),
                         surv=c(surv_A, surv_B),
                         group=c(rep(0, length(times$overall)),
                                 rep(1, length(times$overall))))
  return(plotdata)
}

## Direct Adjustment
#' @export
surv_method_direct <- function(data, variable, ev_time, event, sd,
                               outcome_model, times, verbose=F, ...) {

  surv <- riskRegression::ate(event=outcome_model, treatment=variable,
                              data=data, estimator="Gformula",
                              times=times$overall, se=sd, verbose=verbose, ...)
  plotdata <- data.frame(time=surv$meanRisk$time,
                         surv=1 - surv$meanRisk$estimate,
                         group=surv$meanRisk$treatment)

  if (sd) {
    plotdata$sd <- surv$meanRisk$se * sqrt(nrow(data))
  }

  return(plotdata)
}

## Using propensity score matching
#' @export
surv_method_matching <- function(data, variable, ev_time, event, sd,
                                 treatment_model, ...) {

  if (is.numeric(treatment_model)) {
    ps_score <- treatment_model
  } else {
    ps_score <- predict(treatment_model, newdata=data, type="response")
  }

  rr <- Matching::Match(Tr=data[, variable], X=ps_score, estimand="ATE", ...)
  m_dat <- rbind(data[rr$index.treated,], data[rr$index.control,])

  # estimate survival curve
  form <- paste0("Surv(", ev_time, ", ", event, ") ~ ", variable)
  surv <- survfit(as.formula(form), data=m_dat, se.fit=sd)
  plotdata <- data.frame(time=surv$time,
                         est=surv$surv,
                         var=c(rep(0, surv$strata[1]),
                               rep(1, surv$strata[2])))

  if (sd) {
    plotdata$sd <- surv$std.err * sqrt(nrow(data))
  }

  return(plotdata)
}

## Using Augmented Inverse Probability of Treatment Weighting
#' @export
surv_method_aiptw <- function(data, variable, ev_time, event, sd, times,
                              outcome_model=NULL, treatment_model=NULL,
                              censoring_model=NULL, verbose=F, ...) {

  if (is.null(censoring_model) & is.null(treatment_model) &
      is.null(outcome_model)) {
    stop("At least one of 'treatment_model', 'outcome_model' and ",
         "'censoring_model' needs to be specified, see details.")
  }

  if (is.null(censoring_model)) {
    form <- paste0("Surv(", ev_time, ", ", event, "==0) ~ 1")
    censoring_model <- coxph(as.formula(form), data=data, x=T)
  }
  if (is.null(treatment_model)) {
    form <- paste0(variable, " ~ 1")
    treatment_model <- glm(as.formula(form), data=data, family="binomial")
  }
  if (is.null(outcome_model)) {
    form <- paste0("Surv(", ev_time, ", ", event, ") ~ 1")
    outcome_model <- coxph(as.formula(form), data=data, x=T)
  }

  # estimate AIPTW cumulative incidence
  curve <- riskRegression::ate(event=outcome_model,
                               treatment=treatment_model,
                               censor=censoring_model,
                               data=data,
                               times=times$overall,
                               se=sd,
                               verbose=verbose,
                               estimator="AIPTW,AIPCW")$meanRisk
  # calculate survival estimates
  plotdata <- data.frame(time=c(curve$time[curve$treatment==0],
                                curve$time[curve$treatment==1]),
                         surv=c(1 - curve$estimate[curve$treatment==0],
                                1 - curve$estimate[curve$treatment==1]),
                         group=c(rep(0, length(times$overall)),
                                 rep(1, length(times$overall))))
  if (sd) {
    plotdata$sd <- c(curve$se[curve$treatment==0],
                     curve$se[curve$treatment==1]) * sqrt(nrow(data))
  }

  return(plotdata)
}

## function to get predicted values from any kind of geese object,
## for multiple time points
geese_predictions <- function(geese, Sdata, times, n) {

  current.na.action <- options('na.action')
  options(na.action="na.pass")
  # full model matrix and betas
  mod_mat <- model.matrix(geese$formula, data=Sdata)
  options(na.action=current.na.action[[1]])

  betas <- geese$beta

  apply_betas <- function(x, betas) {
    return(sum(x * betas))
  }

  pred_mat <- matrix(nrow=n, ncol=length(times))
  for (i in 1:length(times)) {
    # take only relevant portion (at time t) of model matrix
    mod_mat_t <- mod_mat[Sdata$vtime==times[i],]
    # apply coefficients
    preds <- apply(X=mod_mat_t, MARGIN=1, FUN=apply_betas, betas=betas)
    pred_mat[,i] <- preds
  }
  colnames(pred_mat) <- paste0("t.", times)
  return(pred_mat)
}

## Using Pseudo Observations and Direct Adjustment
#' @export
surv_method_direct_pseudo <- function(data, variable, ev_time, event, sd,
                                      times, outcome_vars, type_time="factor",
                                      spline_df=10, na.rm=F) {
  # some constants
  len <- length(times$overall)
  n <- nrow(data)
  group <- data[,variable]

  # estimate pseudo observations
  hist_formula <- as.formula(paste("Hist(", ev_time, ", ", event, ") ~ 1"))
  pseudo <- prodlim::jackknife(prodlim(hist_formula, data=data),
                               times=times$overall)
  # create data for geese
  Sdata <- data.frame(yi=c(pseudo),
                      group=rep(group, len),
                      vtime=rep(times$overall, rep(n, len)),
                      id=rep(1:n, len))
  for (col in outcome_vars) {
    Sdata[,col] <- rep(data[,col], len)
  }

  if (type_time=="factor") {
    Sdata$vtime <- as.factor(Sdata$vtime)
    geese_formula <- paste("yi ~ vtime + ", paste(outcome_vars, collapse=" + "),
                           " + group")
  } else if (type_time=="bs") {
    geese_formula <- paste("yi ~ bs(vtime, df=", spline_df, ") + ",
                           paste(outcome_vars, collapse=" + "), " + group")
  } else if (type_time=="ns") {
    geese_formula <- paste("yi ~ ns(vtime, df=", spline_df, ") + ",
                           paste(outcome_vars, collapse=" + "), " + group")
  }

  # call geese
  geese_mod <- geepack::geese(as.formula(geese_formula), scale.fix=TRUE,
                              data=Sdata, family=gaussian, id=id, jack=F,
                              mean.link="cloglog", corstr="independence")

  # for group = 1
  Sdata$group <- 1
  pred_1 <- geese_predictions(geese_mod, Sdata, times=times$overall, n=n)

  m1 <- 1 - exp(-exp(pred_1))
  surv1 <- apply(m1, 2, mean, na.rm=na.rm)

  # for group = 0
  Sdata$group <- 0
  pred_0 <- geese_predictions(geese_mod, Sdata, times=times$overall, n=n)

  m0 <- 1 - exp(-exp(pred_0))
  surv0 <- apply(m0, 2, mean, na.rm=na.rm)

  plotdata <- data.frame(time=c(times$overall, times$overall),
                         surv=c(surv0, surv1),
                         group=c(rep(0, len), rep(1, len)))
  if (sd) {
    pseudo_direct_sd <- function(x, n, na.rm) {
      sqrt(var(x, na.rm=na.rm) / n)
    }

    survsd1 <- apply(m1, 2, pseudo_direct_sd, n=n, na.rm=na.rm)
    survsd0 <- apply(m0, 2, pseudo_direct_sd, n=n, na.rm=na.rm)
    plotdata$sd <- c(survsd0, survsd1)
  }
  return(plotdata)
}

## Using Pseudo Observations Doubly-Robust, Wang (2018)
#' @export
surv_method_aiptw_pseudo <- function(data, variable, ev_time, event, sd,
                                     times, outcome_vars, treatment_model,
                                     type_time="factor", spline_df=10, na.rm=F) {
  # some constants
  len <- length(times$overall)
  n <- nrow(data)
  group <- data[,variable]

  if (is.numeric(treatment_model)) {
    ps_score <- treatment_model
  } else {
    ps_score <- treatment_model$fitted.values
  }

  # estimate pseudo observations
  hist_formula <- as.formula(paste("Hist(", ev_time, ", ", event, ") ~ 1"))
  pseudo <- prodlim::jackknife(prodlim(hist_formula, data=data),
                               times=times$overall)
  # create data for geese
  Sdata <- data.frame(yi=c(pseudo),
                      group=rep(group, len),
                      vtime=rep(times$overall, rep(n, len)),
                      id=rep(1:n, len))
  for (col in outcome_vars) {
    Sdata[,col] <- rep(data[,col], len)
  }

  if (type_time=="factor") {
    Sdata$vtime <- as.factor(Sdata$vtime)
    geese_formula <- paste("yi ~ vtime + ", paste(outcome_vars, collapse=" + "),
                           " + group")
  } else if (type_time=="bs") {
    geese_formula <- paste("yi ~ bs(vtime, df=", spline_df, ") + ",
                           paste(outcome_vars, collapse=" + "), " + group")
  } else if (type_time=="ns") {
    geese_formula <- paste("yi ~ ns(vtime, df=", spline_df, ") + ",
                           paste(outcome_vars, collapse=" + "), " + group")
  }

  # call geese
  geese_mod <- geepack::geese(as.formula(geese_formula), scale.fix=TRUE,
                              data=Sdata, family=gaussian, id=id, jack=F,
                              mean.link="cloglog", corstr="independence")

  # for group = 1
  Sdata$group <- 1
  pred_1 <- geese_predictions(geese_mod, Sdata, times=times$overall, n=n)

  m1 <- 1 - exp(-exp(pred_1))
  dr1 <- (pseudo*group-(group-ps_score)*m1)/ps_score
  surv1 <- apply(dr1, 2, mean, na.rm=na.rm)

  # for group = 0
  Sdata$group <- 0
  pred_0 <- geese_predictions(geese_mod, Sdata, times=times$overall, n=n)

  m0 <- 1 - exp(-exp(pred_0))
  dr0 <- (pseudo*(1-group)+(group-ps_score)*m0)/(1-ps_score)
  surv0 <- apply(dr0, 2, mean, na.rm=na.rm)

  plotdata <- data.frame(time=c(times$overall, times$overall),
                         surv=c(surv0, surv1),
                         group=c(rep(0, len), rep(1, len)))
  if (sd) {
    pseudo_dr_sd <- function(x, n, na.rm) {
      sqrt(var(x, na.rm=na.rm) / n)
    }

    survsd1 <- apply(dr1, 2, pseudo_dr_sd, n=n, na.rm=na.rm)
    survsd0 <- apply(dr0, 2, pseudo_dr_sd, n=n, na.rm=na.rm)
    plotdata$sd <- c(survsd0, survsd1)
  }

  return(plotdata)
}

## Using Empirical Likelihood Estimation
#' @export
surv_method_el <- function(data, variable, ev_time, event, sd,
                           times, treatment_vars, moment="first",
                           standardize=F, ...) {

  el_0 <- adjKMtest::el.est(y=data[, ev_time],
                            delta=data[, event],
                            treat=data[, variable],
                            x=as.matrix(data[, treatment_vars]),
                            treat.select=0,
                            t=times$overall,
                            psix_moment=moment,
                            standardize=standardize,
                            get.sd=sd,
                            ...)
  el_1 <- adjKMtest::el.est(y=data[, ev_time],
                            delta=data[, event],
                            treat=data[, variable],
                            x=as.matrix(data[, treatment_vars]),
                            treat.select=1,
                            t=times$overall,
                            psix_moment=moment,
                            standardize=standardize,
                            get.sd=sd,
                            ...)

  plotdata <- data.frame(time=c(times$overall, times$overall),
                         surv=c(el_0$St, el_1$St),
                         group=c(rep(0, length(el_0$St)),
                                 rep(1, length(el_0$St))))
  if (sd) {
    plotdata$sd <- c(el_0$sd, el_1$sd)
  }

  return(plotdata)
}

## Targeted Maximum Likelihood Estimation
#' @export
surv_method_tmle <- function(data, variable, ev_time, event, sd,
                             times, t_max=NULL, adjustVars=NULL,
                             SL.ftime=NULL, SL.ctime=NULL, SL.trt=NULL,
                             glm.ftime=NULL, glm.ctime=NULL, glm.trt=NULL,
                             ...) {
  if (is.null(t_max)) {
    t_max <- max(times$overall)
  }
  if (is.null(adjustVars)) {
    all_covars <- colnames(data)
    all_covars <- all_covars[!all_covars %in% c(variable, ev_time, event)]
    adjustVars <- data[, all_covars]
  }

  # TMLE fit
  fit <- survtmle::survtmle(
    ftime=data[, ev_time],
    ftype=data[, event],
    trt=data[, variable],
    t0=t_max,
    adjustVars=adjustVars,
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
    survtmle::timepoints(fit, times=times$overall)
  }, warning=function(w) {
    if (startsWith(conditionMessage(w), "Using formula(x) is deprecated"))
      invokeRestart("muffleWarning")
  })
  len_groups <- length(unique(data[, variable]))

  # extract estimates of survival
  names_groups <- unique(
    lapply(lapply(tpfit, FUN=`[[`, "est"), FUN=rownames)
  )[[1]]

  est_only <- t(matrix(unlist(
    lapply(tpfit, FUN=`[[`, "est")
  ), ncol=len_groups, byrow=TRUE))
  est_only <- as.data.frame(est_only)
  rownames(est_only) <- names_groups

  s_0 <- 1 - as.numeric(est_only[1, ])
  s_1 <- 1 - as.numeric(est_only[2, ])

  # put together
  plotdata <- data.frame(time=rep(times$overall, 2),
                         surv=c(s_0, s_1),
                         group=c(rep(0, length(times$overall)),
                                 rep(1, length(times$overall))))

  return(plotdata)
}

## One-Step Targeted Maximum Likelihood Estimation
#' @export
surv_method_ostmle <- function(data, variable, ev_time, event, sd,
                               times, t_max=NULL, adjust_vars=NULL,
                               SL.ftime=NULL, SL.ctime=NULL, SL.trt=NULL,
                               epsilon=1, max_num_iteration=100,
                               psi_moss_method="l2", tmle_tolerance=NULL,
                               gtol=1e-3) {
  if (is.null(t_max)) {
    t_max <- max(times$overall)
  }
  if (is.null(adjust_vars)) {
    all_covars <- colnames(data)
    all_covars <- all_covars[!all_covars %in% c(variable, ev_time, event)]
    adjust_vars <- data[, all_covars]
  }

  # time point grid
  k_grid <- 1:max(data[, ev_time])

  # create initial fit object
  sl_fit <- MOSS::initial_sl_fit(
    T_tilde=data[, ev_time],
    Delta=data[, event],
    A=data[, variable],
    W=adjust_vars,
    t_max=t_max,
    sl_treatment=SL.trt,
    sl_censoring=SL.ctime,
    sl_failure=SL.ftime,
    gtol=gtol
  )
  invisible(sl_fit$density_failure_1$hazard_to_survival())
  invisible(sl_fit$density_failure_0$hazard_to_survival())

  sl_fit$density_failure_1$t <- 1:t_max
  sl_fit$density_failure_0$t <- 1:t_max

  moss_hazard_fit <- MOSS::MOSS_hazard$new(
    A=data[, variable],
    T_tilde=data[, ev_time],
    Delta=data[, event],
    density_failure=sl_fit$density_failure_1,
    density_censor=sl_fit$density_censor_1,
    g1W=sl_fit$g1W,
    A_intervene=1,
    k_grid=k_grid
  )
  psi_moss_hazard_1 <- moss_hazard_fit$iterate_onestep(
    epsilon=epsilon,
    max_num_interation=max_num_iteration,
    tmle_tolerance=tmle_tolerance,
    verbose=FALSE,
    method=psi_moss_method
  )

  moss_hazard_fit <- MOSS::MOSS_hazard$new(
    A=data[, variable],
    T_tilde=data[, ev_time],
    Delta=data[, event],
    density_failure=sl_fit$density_failure_0,
    density_censor=sl_fit$density_censor_0,
    g1W=sl_fit$g1W,
    A_intervene=0,
    k_grid=k_grid
  )
  psi_moss_hazard_0 <- moss_hazard_fit$iterate_onestep(
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
  return(plotdata)
}
