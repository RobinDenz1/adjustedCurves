
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
  times <- sort(unique(adjsurv$adj$time))
  n.risk <- get_risk_table(times=times,
                           data=adjsurv$data,
                           ev_time=adjsurv$ev_time,
                           variable=adjsurv$variable,
                           event=adjsurv$event,
                           type="n_at_risk",
                           weights=weights,
                           digits=Inf)
  colnames(n.risk)[colnames(n.risk)=="group"] <- "strata"
  colnames(n.risk)[colnames(n.risk)=="est"] <- "n.risk"

  n.risk$strata <- factor(n.risk$strata, levels=levels(df$strata))
  df <- merge(df, n.risk, by=c("time", "strata"), all.x=TRUE, all.y=FALSE)

  n.event <- get_risk_table(times=times,
                            data=adjsurv$data,
                            ev_time=adjsurv$ev_time,
                            variable=adjsurv$variable,
                            event=adjsurv$event,
                            type="raw_events",
                            weights=weights,
                            digits=Inf)
  colnames(n.event)[colnames(n.event)=="group"] <- "strata"
  colnames(n.event)[colnames(n.event)=="est"] <- "n.event"

  n.event$strata <- factor(n.event$strata, levels=levels(df$strata))
  df <- merge(df, n.event, by=c("time", "strata"))

  df <- df[order(df$strata, df$time), ]

  return(df)
}
