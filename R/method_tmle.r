
## redefine 'timepoints' from survtmle to fix a bug in there
survtmle.timepoints <- function(object, times, returnModels=FALSE,
                                SL.trt, SL.ctime, SL.ftime,
                                glm.trt, glm.ctime, glm.ftime) {

  callList <- as.list(object$call)[-1]
  cglm <- any(class(object$ctimeMod) %in% c("glm", "speedglm")) |
    any(class(object$ctimeMod) == "noCens")

  tglm <- any(class(object$trtMod) %in% c("glm", "speedglm"))
  ftglm <- ifelse(callList$method == "hazard",
                  any(class(object$ftimeMod[[1]]) %in% c(
                    "glm",
                    "speedglm"
                  )), FALSE
  )

  myOpts <- c(
    "t0", "returnModels",
    ifelse(cglm, "glm.ctime", "SL.ctime"),
    ifelse(tglm, "glm.trt", "SL.trt")
  )
  if (callList$method == "hazard") {
    myOpts <- c(myOpts, ifelse(ftglm, "glm.ftime", "SL.ftime"))
  }
  funOpts <- callList[-which(names(callList) %in% myOpts)]

  funOpts$returnModels <- returnModels
  # used glm for censoring?
  if (cglm) {
    funOpts$glm.ctime <- object$ctimeMod
    funOpts$SL.ctime <- NULL
  } else {
    funOpts$SL.ctime <- object$ctimeMod
  }
  # used glm for trt?
  if (tglm) {
    funOpts$glm.trt <- object$trtMod
  } else {
    funOpts$SL.trt <- object$trtMod
  }
  # used glm for ftime
  if (ftglm & callList$method == "hazard") {
    funOpts$glm.ftime <- object$ftimeMod
  } else if (!ftglm & callList$method == "hazard") {
    funOpts$SL.ftime <- object$ftimeMod
  }
  # NOTE: this is the bug-fix
  # add in failure times, types, trt, and adjust
  funOpts$ftime <- object$ftime
  funOpts$ftype <- object$ftype
  funOpts$trt <- object$trt
  funOpts$adjustVars <- object$adjustVars
  funOpts$ftypeOfInterest <- object$ftypeOfInterest

  outList <- vector(mode = "list", length = length(times))
  ct <- 0
  for (i in times) {
    ct <- ct + 1
    funOpts$t0 <- i
    if (all(object$ftime[object$ftype > 0] > i)) {
      outList[[ct]] <- list(
        est = rep(0, length(object$est)),
        var = matrix(
          NA,
          nrow = length(object$est),
          ncol = length(object$est)
        )
      )
    } else {
      if (i != object$t0) {
        outList[[ct]] <- do.call("survtmle", args = funOpts)
      } else {
        outList[[ct]] <- object
      }
    }
  }
  names(outList) <- paste0("t", times)
  class(outList) <- "tp.survtmle"
  return(outList)
}

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
