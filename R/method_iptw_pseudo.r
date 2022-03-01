
## generalized function that does both
method_iptw_pseudo <- function(data, variable, ev_time, event, cause, conf_int,
                               conf_level=0.95, times, treatment_model,
                               weight_method="ps", stabilize=FALSE,
                               trim=FALSE, se_method="cochrane",
                               censoring_vars=NULL, ipcw_method="binder",
                               mode, ...) {

  levs <- levels(data[, variable])

  # get weights
  if (is.numeric(treatment_model)) {
    weights <- treatment_model
    weights <- trim_weights(weights, trim)
    if (stabilize) {
      weights <- stabilize_weights(weights, data, variable, levs)
    }
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
                             ipcw.method=ipcw_method,
                             cause=cause)

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

  if (mode=="cif") {
    colnames(plotdata)[colnames(plotdata)=="surv"] <- "cif"
  }

  output <- list(plotdata=plotdata,
                 pseudo_values=pseudo,
                 weights=weights)

  if (mode=="cif") {
    class(output) <- "adjustedcif.method"
  } else {
    class(output) <- "adjustedsurv.method"
  }

  return(output)
}

## Using Pseudo-Observations and IPTW for Survival Curves
#' @export
surv_iptw_pseudo <- function(data, variable, ev_time, event, conf_int,
                             conf_level=0.95, times, treatment_model,
                             weight_method="ps", stabilize=FALSE,
                             trim=FALSE, se_method="cochrane",
                             censoring_vars=NULL, ipcw_method="binder", ...) {

  out <- method_iptw_pseudo(data=data, variable=variable, ev_time=ev_time,
                            event=event, conf_int=conf_int,
                            conf_level=conf_level, times=times,
                            treatment_model=treatment_model,
                            weight_method=weight_method, stabilize=stabilize,
                            trim=trim, se_method=se_method,
                            censoring_vars=censoring_vars,
                            ipcw_method=ipcw_method, cause=1,
                            mode="surv", ...)
  return(out)
}

## Using Pseudo-Observations and IPTW for CIFs
#' @export
cif_iptw_pseudo <- function(data, variable, ev_time, event, cause,
                            conf_int, conf_level=0.95, times,
                            treatment_model, weight_method="ps",
                            stabilize=FALSE, trim=FALSE,
                            se_method="cochrane", ...) {

  out <- method_iptw_pseudo(data=data, variable=variable, ev_time=ev_time,
                            event=event, cause=cause, conf_int=conf_int,
                            conf_level=conf_level, times=times,
                            treatment_model=treatment_model,
                            weight_method=weight_method, stabilize=stabilize,
                            trim=trim, se_method=se_method,
                            censoring_vars=NULL, ipcw_method="binder",
                            mode="cif", ...)
  return(out)
}
