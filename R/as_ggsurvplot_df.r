
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

  if (adjsurv$method=="iptw_km") {
    df$n.risk <- adjsurv$n_at_risk$n_at_risk
    df$n.event <- adjsurv$n_at_risk$n_events
  } else if (adjsurv$method=="km") {
    df$n.risk <- adjsurv$survfit_object$n.risk
    df$n.event <- adjsurv$survfit_object$n.event
  }

  return(df)
}
