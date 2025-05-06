
## Using Propensity Score Matching
#' @export
surv_matching <- function(data, variable, ev_time, event, conf_int=FALSE,
                          conf_level=0.95, times, treatment_model,
                          gtol=0.001, ...) {

  # turn it into numeric
  levs <- levels(data[, variable])
  data[, variable] <- ifelse(data[,variable]==levs[1], 0, 1)

  if (is.numeric(treatment_model)) {
    ps_score <- treatment_model
  } else {
    ps_score <- stats::predict.glm(treatment_model, newdata=data,
                                   type="response")
  }

  # trim extreme propensity score according to gtol
  ps_score[ps_score < gtol] <- gtol
  ps_score[ps_score > (1 - gtol)] <- 1 - gtol

  # perform matching
  rr <- Matching::Match(Tr=data[, variable], X=ps_score, estimand="ATE", ...)
  m_dat <- rbind(data[rr$index.treated, ], data[rr$index.control, ])

  weights <- rr$weights
  m_dat$match_weights <- c(weights, weights)

  # estimate survival curve
  form <- paste0("survival::Surv(", ev_time, ", ", event, ") ~ ", variable)
  surv <- survival::survfit(stats::as.formula(form), data=m_dat,
                            se.fit=conf_int, conf.int=conf_level,
                            weights=m_dat$match_weights,
                            robust=TRUE)
  plotdata <- data.frame(time=surv$time,
                         surv=surv$surv,
                         group=c(rep(levs[1], surv$strata[1]),
                                 rep(levs[2], surv$strata[2])))

  if (!is.null(times)) {
    plotdata <- specific_times(plotdata, times)
  }
  plotdata$group <- factor(plotdata$group, levels=levs)

  output <- list(plotdata=plotdata,
                 match_object=rr,
                 survfit_object=surv)
  class(output) <- "adjustedsurv.method"

  return(output)
}

## Matching
#' @export
cif_matching <- function(data, variable, ev_time, event, cause, conf_int,
                         conf_level=0.95, times, treatment_model,
                         gtol=0.001, ...) {

  # turn it into numeric
  levs <- levels(data[, variable])
  data[, variable] <- ifelse(data[, variable]==levs[1], 0, 1)

  if (is.numeric(treatment_model)) {
    ps_score <- treatment_model
  } else {
    ps_score <- stats::predict.glm(treatment_model, newdata=data,
                                   type="response")
  }

  # trim extreme propensity score according to gtol
  ps_score[ps_score < gtol] <- gtol
  ps_score[ps_score > (1 - gtol)] <- 1 - gtol

  # perform matching
  rr <- Matching::Match(Tr=data[, variable], X=ps_score, estimand="ATE", ...)
  m_dat <- rbind(data[rr$index.treated, ], data[rr$index.control, ])

  # estimate cif
  plotdata <- cif_aalen_johansen(data=m_dat,
                                 variable=variable,
                                 ev_time=ev_time,
                                 event=event,
                                 cause=cause,
                                 conf_int=FALSE,
                                 conf_level=conf_level)$plotdata

  if (!is.null(times)) {
    plotdata$se <- NULL
    plotdata <- specific_times(plotdata, times, est="cif")
  }

  # get factor levels back
  plotdata$group[plotdata$group==0] <- levs[1]
  plotdata$group[plotdata$group==1] <- levs[2]
  plotdata$group <- factor(plotdata$group, levels=levs)

  output <- list(plotdata=plotdata,
                 match_object=rr)
  class(output) <- "adjustedcif.method"

  return(output)
}
