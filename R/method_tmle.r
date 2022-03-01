
## Targeted Maximum Likelihood Estimation
#' @export
cif_tmle <- function(data, variable, ev_time, event, cause, conf_int,
                     conf_level=0.95, times, adjust_vars=NULL,
                     SL.ftime=NULL, SL.ctime=NULL, SL.trt=NULL,
                     glm.ftime=NULL, glm.ctime=NULL, glm.trt=NULL,
                     ...) {
  # if it's a factor, turn it into numeric
  if (is.factor(data[, variable])) {
    levs <- levels(data[, variable])
    data[, variable] <- ifelse(data[, variable]==levs[1], 0, 1)
  } else {
    levs <- unique(data[, variable])
  }

  # gather needed data
  if (is.null(adjust_vars)) {
    all_covars <- colnames(data)
    all_covars <- all_covars[!all_covars %in% c(variable, ev_time, event)]
    adjust_vars <- data[, all_covars]
  } else {
    adjust_vars <- data[, adjust_vars]
  }

  # TMLE fit
  fit <- survtmle::survtmle(
    ftime=data[, ev_time],
    ftype=data[, event],
    trt=data[, variable],
    t0=max(data[, ev_time]),
    adjustVars=adjust_vars,
    ftypeOfInterest=cause,
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
    survtmle.timepoints(fit, times=times, SL.trt=SL.trt, SL.ctime=SL.ctime,
                        SL.ftime=SL.ftime, glm.trt=glm.trt,
                        glm.ctime=glm.ctime, glm.ftime=glm.ftime)
  }, warning=function(w) {
    if (startsWith(conditionMessage(w), "Using formula(x) is deprecated"))
      invokeRestart("muffleWarning")
  })

  cif_0 <- unlist(lapply(tpfit, function(x) {x$est[1]}))
  cif_1 <- unlist(lapply(tpfit, function(x) {x$est[2]}))

  # put together
  plotdata <- data.frame(time=rep(times, 2),
                         cif=c(cif_0, cif_1),
                         group=c(rep(levs[1], length(times)),
                                 rep(levs[2], length(times))))

  if (conf_int) {

    var_0 <- unlist(lapply(tpfit, function(x) {x$var[1]}))
    var_1 <- unlist(lapply(tpfit, function(x) {x$var[2]}))

    plotdata$se <- sqrt(c(var_0, var_1))

    confint.tp.survtmle <- utils::getFromNamespace("confint.tp.survtmle",
                                                   "survtmle")
    survtmle_ci <- confint.tp.survtmle(tpfit, level=conf_level)

    plotdata$ci_lower <- c(survtmle_ci$`0 1`[, 1], survtmle_ci$`1 1`[, 1])
    plotdata$ci_upper <- c(survtmle_ci$`0 1`[, 2], survtmle_ci$`1 1`[, 2])

  }

  output <- list(plotdata=plotdata,
                 survtmle_object=fit,
                 survtmle.timepoints_object=tpfit)
  class(output) <- "adjustedcif.method"

  return(output)
}

## Targeted Maximum Likelihood Estimation
#' @export
surv_tmle <- function(data, variable, ev_time, event, conf_int,
                      conf_level=0.95, times, adjust_vars=NULL,
                      SL.ftime=NULL, SL.ctime=NULL, SL.trt=NULL,
                      glm.ftime=NULL, glm.ctime=NULL, glm.trt=NULL,
                      ...) {

  out <- cif_tmle(data=data, variable=variable, ev_time=ev_time,
                  event=event, conf_int=conf_int, conf_level=conf_level,
                  times=times, adjust_vars=adjust_vars,
                  SL.ftime=SL.ftime, SL.ctime=SL.ctime, SL.trt=SL.trt,
                  glm.ftime=glm.ftime, glm.ctime=glm.ctime,
                  glm.trt=glm.trt, cause=1, ...)

  plotdata <- out$plotdata
  colnames(plotdata)[colnames(plotdata)=="cif"] <- "surv"

  plotdata$surv <- 1 - plotdata$surv

  if ("ci_lower" %in% colnames(plotdata)) {
    plotdata$ci_lower <- 1 - plotdata$ci_lower
    plotdata$ci_upper <- 1 - plotdata$ci_upper
  }

  out$plotdata <- plotdata
  class(out) <- "adjustedsurv.method"

  return(out)
}
