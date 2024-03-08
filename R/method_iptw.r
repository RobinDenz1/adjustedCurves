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

## IPTW Kaplan-Meier estimate
#' @export
surv_iptw_km <- function(data, variable, ev_time, event, conf_int,
                         conf_level=0.95, times=NULL, treatment_model,
                         weight_method="ps", stabilize=FALSE,
                         trim=FALSE, trim_quantiles=FALSE, ...) {

  levs <- levels(data[, variable])

  # get weights
  if (is.numeric(treatment_model)) {
    weights <- treatment_model
    weights <- trim_weights(weights=weights, trim=trim)
    weights <- trim_weights_quantiles(weights=weights, trim_q=trim_quantiles)
    if (stabilize) {
      weights <- stabilize_weights(weights, data, variable, levs)
    }
  } else {
    weights <- get_iptw_weights(data=data, treatment_model=treatment_model,
                                weight_method=weight_method,
                                variable=variable, stabilize=stabilize,
                                trim=trim, trim_q=trim_quantiles, ...)
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

    plotdata_group <- data.frame(time=tj, surv=st, group=levs[i],
                                 n_risk=yj, n_event=dj)

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

  # extract weighted n at risk
  n_at_risk <- data.frame(time=plotdata$time,
                          group=plotdata$group,
                          n_at_risk=plotdata$n_risk,
                          n_events=plotdata$n_event)
  plotdata$n_risk <- NULL
  plotdata$n_event <- NULL

  if (!is.null(times)) {
    plotdata <- specific_times(plotdata, times)
  }

  output <- list(plotdata=plotdata,
                 weights=weights,
                 n_at_risk=n_at_risk)
  class(output) <- "adjustedsurv.method"

  return(output)
}

## IPTW with univariate cox-model
# NOTE:
# - using the G-Formula directly on the iptw cox model does not work,
#   probably because all predict() methods ignore the weights when
#   calculating the baseline hazard. Only this version is actually unbiased
#' @export
surv_iptw_cox <- function(data, variable, ev_time, event, conf_int=FALSE,
                          conf_level=0.95, times=NULL, treatment_model,
                          weight_method="ps", stabilize=FALSE,
                          trim=FALSE, trim_quantiles=FALSE, ...) {

  levs <- levels(data[, variable])

  # get weights
  if (is.numeric(treatment_model)) {
    weights <- treatment_model
    weights <- trim_weights(weights=weights, trim=trim)
    weights <- trim_weights_quantiles(weights=weights, trim_q=trim_quantiles)
    if (stabilize) {
      weights <- stabilize_weights(weights, data, variable, levs)
    }
  } else {
    weights <- get_iptw_weights(data=data, treatment_model=treatment_model,
                                weight_method=weight_method,
                                variable=variable, stabilize=stabilize,
                                trim=trim, trim_q=trim_quantiles, ...)
  }

  # univariate, weighted cox model
  form <- paste0("survival::Surv(", ev_time, ", ", event, ") ~ strata(",
                 variable, ")")
  model <- survival::coxph(stats::as.formula(form), weights=weights, data=data)
  surv <- survival::survfit(model, se.fit=FALSE, conf.int=conf_level)

  plotdata <- data.frame(time=surv$time,
                         surv=surv$surv)

  # get grouping variable
  group <- c()
  for (strat in names(surv$strata)) {
    group <- c(group, rep(strat, surv$strata[strat]))
  }
  group <- gsub(paste0(variable, "="), "", group)
  plotdata$group <- group

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
