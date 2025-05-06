
## simple Kaplan-Meier estimate
#' @export
surv_km <- function(data, variable, ev_time, event, conf_int,
                    conf_level=0.95, times=NULL, conf_type="log") {

  form <- paste0("survival::Surv(", ev_time, ", ", event, ") ~ ", variable)

  surv <- survival::survfit.formula(stats::as.formula(form), data=data,
                                    se.fit=conf_int, conf.int=conf_level,
                                    conf.type=conf_type)
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

  if (!is.null(times)) {
    plotdata <- specific_times(plotdata, times)
  }

  output <- list(plotdata=plotdata,
                 survfit_object=surv)
  class(output) <- "adjustedsurv.method"

  return(output)
}

## Aalen-Johansen estimator
#' @export
cif_aalen_johansen <- function(data, variable, ev_time, event, cause,
                               conf_int, conf_level=0.95, times=NULL, ...) {

  cif <- cmprsk::cuminc(ftime=data[, ev_time],
                        fstatus=data[, event],
                        group=data[, variable],
                        ...)

  levs <- unique(data[, variable])
  cif_names <- paste(levs, cause)
  plotdata <- vector(mode="list", length=length(cif_names))
  for (i in seq_len(length(cif_names))) {
    plotdata[[i]] <- data.frame(time=cif[[cif_names[i]]]$time,
                                cif=cif[[cif_names[i]]]$est,
                                group=levs[i],
                                se=cif[[cif_names[i]]]$var)
  }
  plotdata <- as.data.frame(dplyr::bind_rows(plotdata))

  if (conf_int) {
    plotdata$se <- sqrt(plotdata$se)
    cif_cis <- confint_surv(surv=plotdata$cif, se=plotdata$se,
                            conf_level=conf_level, conf_type="plain")
    plotdata$ci_lower <- cif_cis$left
    plotdata$ci_upper <- cif_cis$right
  } else {
    plotdata$se <- NULL
  }

  # remove weird structure in cmprsk::cuminc call
  ids <- seq(2, nrow(plotdata), 2)
  ids <- ids[1:(length(ids-1))]
  plotdata <- plotdata[ids, ]

  if (!is.null(times)) {
    plotdata <- specific_times(plotdata, times, est="cif")
  }

  output <- list(plotdata=plotdata,
                 cuminc_object=cif)
  class(output) <- "adjustedcif.method"

  return(output)
}
