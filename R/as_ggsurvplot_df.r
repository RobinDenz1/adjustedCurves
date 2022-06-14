
## Get a data.frame from adjustedsurv objects that can be
## used directly to call the ggsurvplot_df() function from the
## survminer package
#' @export
as_ggsurvplot_df <- function(adjsurv) {

  if (!inherits(adjsurv, "adjustedsurv")) {
    stop("This function can only be used with 'adjustedsurv' objects,",
         " created using the 'adjustedsurv' function.")
  }

  df <- data.frame(time=adjsurv$adjsurv$time,
                   surv=adjsurv$adjsurv$surv,
                   strata=adjsurv$adjsurv$group)

  if ("se" %in% colnames(adjsurv)) {
    df$std.err <- adjsurv$adjsurv$se
    df$upper <- adjsurv$adjsurv$ci_upper
    df$lower <- adjsurv$adjsurv$ci_lower
  }

  if (adjsurv$method=="iptw_km") {
    df$n.risk <- adjsurv$n_at_risk$n_at_risk
    df$n.event <- adjsurv$n_at_risk$n_events
  }

  return(df)
}
