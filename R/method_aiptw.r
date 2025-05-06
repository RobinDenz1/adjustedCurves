
## Using Augmented Inverse Probability of Treatment Weighting for CIFs
#' @export
cif_aiptw <- function(data, variable, ev_time, event, cause, conf_int,
                      conf_level=0.95, times, outcome_model,
                      treatment_model, censoring_model=NULL,
                      verbose=FALSE, ...) {

  # defaults for input models
  if (is.null(censoring_model)) {
    form <- paste0("survival::Surv(", ev_time, ", ", event, "==0) ~ 1")
    censoring_model <- survival::coxph(stats::as.formula(form), data=data,
                                       x=TRUE, y=TRUE)
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
                               cause=cause,
                               ...)
  # calculate CIF estimates
  plotdata <- data.frame(time=curve$meanRisk$time,
                         cif=curve$meanRisk$estimate,
                         group=curve$meanRisk$treatment)
  if (conf_int) {
    plotdata$se <- curve$meanRisk$se

    confint.ate <- utils::getFromNamespace("confint.ate", "riskRegression")

    cis <- confint.ate(curve, level=conf_level, ci=TRUE)$meanRisk
    plotdata$ci_lower <- cis$lower
    plotdata$ci_upper <- cis$upper
  }

  output <- list(plotdata=plotdata,
                 ate_object=curve)
  class(output) <- "adjustedcif.method"

  return(output)
}

## Using Augmented Inverse Probability of Treatment Weighting for
## Survival Curves
#' @export
surv_aiptw <- function(data, variable, ev_time, event, conf_int,
                       conf_level=0.95, times, outcome_model,
                       treatment_model, censoring_model=NULL,
                       verbose=FALSE, ...) {

  out <- cif_aiptw(data=data, variable=variable, ev_time=ev_time,
                   event=event, conf_int=conf_int, conf_level=conf_level,
                   times=times, outcome_model=outcome_model,
                   treatment_model=treatment_model,
                   censoring_model=censoring_model, verbose=verbose,
                   cause=1, ...)

  plotdata <- out$plotdata
  colnames(plotdata)[colnames(plotdata)=="cif"] <- "surv"

  plotdata$surv <- 1 - plotdata$surv

  if ("ci_lower" %in% colnames(plotdata)) {
    upper <- 1 - plotdata$ci_lower
    lower <- 1 - plotdata$ci_upper
    plotdata$ci_lower <- lower
    plotdata$ci_upper <- upper
  }

  out$plotdata <- plotdata
  class(out) <- "adjustedsurv.method"

  return(out)
}
