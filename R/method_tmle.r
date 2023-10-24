
## targeted maximum likelihood estimation for continuously distributed
## time-to-event data (both with and without competing events)
## as implemented in the concrete R-package
method_tmle <- function(data, variable, ev_time, event, cause, conf_int,
                        conf_level=0.95, times, outcome_model,
                        treatment_model, censoring_model=NULL,
                        cv_args=list(V=5), max_update_iter=500,
                        one_step_eps=0.1,
                        min_nuisance=5/sqrt(nrow(data))/log(nrow(data)),
                        verbose=FALSE, return_models=TRUE) {

  # coerce variable to numeric
  levs <- levels(data[, variable])
  data[, variable] <- ifelse(data[, variable]==levs[1], 0, 1)

  # random censoring by default
  if (is.null(censoring_model)) {
    censoring_model <- list(Surv(time, status==0) ~ 1)
  }

  # construct model list
  model_list <- list("trt" = treatment_model,
                     "0" = censoring_model,
                     "1" = outcome_model)
  names(model_list) <- c(variable, "0", toString(cause))

  # coerce to data.table
  data.table::setDT(data)

  # format arguments
  args <- concrete::formatArguments(DataTable=data,
                                    EventTime=ev_time,
                                    EventType=event,
                                    Treatment=variable,
                                    ID=NULL,
                                    TargetTime=times,
                                    TargetEvent=cause,
                                    Intervention=concrete::makeITT(),
                                    Model=model_list,
                                    CVArg=cv_args,
                                    MaxUpdateIter=max_update_iter,
                                    OneStepEps=one_step_eps,
                                    MinNuisance=min_nuisance,
                                    Verbose=verbose,
                                    GComp=FALSE,
                                    ReturnModels=return_models,
                                    ConcreteArgs=NULL,
                                    RenameCovs=TRUE)

  # do the estimation
  if (verbose) {
    concrete_obj <- concrete::doConcrete(args)
  } else {
    concrete_obj <- suppressMessages(
      concrete::doConcrete(args)
    )
  }

  # get relevant output
  concrete_est <- concrete::getOutput(ConcreteEst=concrete_obj,
                                      Estimand="Risk",
                                      Intervention=seq_along(concrete_obj),
                                      GComp=FALSE,
                                      Simultaneous=FALSE,
                                      Signif=1 - conf_level)

  # extract only neccessary stuff
  concrete_est <- as.data.frame(concrete_est)

  concrete_est <- concrete_est[concrete_est$Event==cause, ]
  concrete_est$Estimand <- NULL
  concrete_est$Estimator <- NULL
  concrete_est$Event <- NULL

  concrete_est$Intervention <- substr(concrete_est$Intervention, 3, 99999)

  # construct output data.frame
  plotdata <- data.frame(time=concrete_est$Time,
                         est=concrete_est$`Pt Est`,
                         group=concrete_est$Intervention)

  # coerce variable back to factor
  plotdata$group <- factor(ifelse(plotdata$group==0, levs[1], levs[2]),
                           levels=levs)

  # whether to include CI or not
  if (conf_int) {
    plotdata$se <- concrete_est$se
    plotdata$ci_lower <- concrete_est$`CI Low`
    plotdata$ci_upper <- concrete_est$`CI Hi`
  }

  # sort it
  plotdata <- plotdata[order(plotdata$group, plotdata$time), ]

  # construct complete output object
  output <- list(plotdata=plotdata,
                 concrete_object=concrete_obj)

  return(output)
}

## TMLE wrapper for simple survival data
#' @export
surv_tmle <- function(data, variable, ev_time, event, conf_int,
                      conf_level=0.95, times, outcome_model,
                      treatment_model, censoring_model=NULL,
                      cv_args=list(V=5), max_update_iter=500,
                      one_step_eps=0.1,
                      min_nuisance=5/sqrt(nrow(data))/log(nrow(data)),
                      verbose=FALSE, return_models=FALSE) {

  out <- method_tmle(data=data, variable=variable, ev_time=ev_time,
                     event=event, cause=1, conf_int=conf_int,
                     conf_level=conf_level, times=times,
                     outcome_model=outcome_model,
                     treatment_model=treatment_model,
                     censoring_model=censoring_model,
                     cv_args=cv_args, max_update_iter=max_update_iter,
                     one_step_eps=one_step_eps,
                     min_nuisance=min_nuisance,
                     verbose=verbose, return_models=return_models)

  # change output data accordingly
  colnames(out$plotdata)[colnames(out$plotdata)=="est"] <- "surv"

  out$plotdata$surv <- 1 - out$plotdata$surv

  if (conf_int) {
    ci_lower <- 1 - out$plotdata$ci_upper
    ci_upper <- 1 - out$plotdata$ci_lower

    out$plotdata$ci_lower <- ci_lower
    out$plotdata$ci_upper <- ci_upper
  }

  class(out) <- "adjustedsurv.method"

  return(out)
}

## TMLE wrapper for competing events data
#' @export
cif_tmle <- function(data, variable, ev_time, event, cause, conf_int,
                     conf_level=0.95, times, outcome_model,
                     treatment_model, censoring_model=NULL,
                     cv_args=list(V=5), max_update_iter=500,
                     one_step_eps=0.1,
                     min_nuisance=5/sqrt(nrow(data))/log(nrow(data)),
                     verbose=FALSE, return_models=FALSE) {

  out <- method_tmle(data=data, variable=variable, ev_time=ev_time,
                     event=event, cause=cause, conf_int=conf_int,
                     conf_level=conf_level, times=times,
                     outcome_model=outcome_model,
                     treatment_model=treatment_model,
                     censoring_model=censoring_model,
                     cv_args=cv_args, max_update_iter=max_update_iter,
                     one_step_eps=one_step_eps,
                     min_nuisance=min_nuisance,
                     verbose=verbose, return_models=return_models)

  colnames(out$plotdata)[colnames(out$plotdata)=="est"] <- "cif"
  class(out) <- "adjustedcif.method"

  return(out)
}
