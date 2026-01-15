
## Get a data.frame from adjustedsurv objects that can be
## used directly to call the ggsurvplot_df() function from the
## survminer package
#' @export
as_ggsurvplot_df <- function(adjsurv) {

  if (!inherits(adjsurv, "adjustedsurv")) {
    stop("This function can only be used with 'adjustedsurv' objects,",
         " created using the 'adjustedsurv' function.")
  }

  df <- data.frame(time=adjsurv$adj$time,
                   surv=adjsurv$adj$surv,
                   strata=adjsurv$adj$group)

  if ("se" %in% colnames(adjsurv$adj)) {
    df$std.err <- adjsurv$adj$se
    df$upper <- adjsurv$adj$ci_upper
    df$lower <- adjsurv$adj$ci_lower
  }

  if (!is.null(adjsurv$mids_analyses)) {

    # set correct weights if specified
    if (!is.null(adjsurv$weights) && is.null(adjsurv$mids_analyses)) {
      weights <- adjsurv$weights
    } else if (!is.null(adjsurv$mids_analyses) &&
               !is.null(adjsurv$mids_analyses[[1]]$weights)) {
      weights <- lapply(adjsurv$mids_analyses, FUN=function(d){d$weights})
    } else {
      weights <- NULL
    }

    # calculate pooled risk table
    n.risk <- get_risk_table(times=df$time,
                             data=adjsurv$data,
                             ev_time=adjsurv$ev_time,
                             variable=adjsurv$variable,
                             event=adjsurv$event,
                             type="n_at_risk",
                             weights=weights,
                             digits=Inf)
    colnames(n.risk) <- c("time", "strata", "n.risk")

    n.risk$n.event <- get_risk_table(times=df$time,
                                     data=adjsurv$data,
                                     ev_time=adjsurv$ev_time,
                                     variable=adjsurv$variable,
                                     event=adjsurv$event,
                                     type="n_events",
                                     weights=weights,
                                     digits=Inf)$est
    df <- merge(df, n.risk, by=c("time", "strata"))

  } else {
    if (adjsurv$method=="iptw_km") {
      df$n.risk <- adjsurv$n_at_risk$n_at_risk
      df$n.event <- adjsurv$n_at_risk$n_events
    } else if (adjsurv$method=="km") {
      df$n.risk <- adjsurv$survfit_object$n.risk
      df$n.event <- adjsurv$survfit_object$n.event
    }
  }

  return(df)
}
