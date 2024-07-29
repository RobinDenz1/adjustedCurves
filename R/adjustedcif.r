# Copyright (C) 2021  Robin Denz
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Second main function of the package. Is basically a wrapper around
## all cif_method functions, offering additional high level stuff
#' @importFrom dplyr %>%
#' @importFrom doRNG %dorng%
#' @importFrom foreach %dopar%
#' @export
adjustedcif <- function(data, variable, ev_time, event, cause, method,
                        conf_int=FALSE, conf_level=0.95, times=NULL,
                        bootstrap=FALSE, n_boot=500, n_cores=1,
                        na.action=options()$na.action,
                        clean_data=TRUE, iso_reg=FALSE,
                        force_bounds=FALSE, mi_extrapolation=FALSE,
                        ...) {

  var_w <- var_b <- B <- var_t <- cif_est <- NULL

  # use data.frame methods only, no tibbles etc.
  if (inherits(data, "data.frame")) {
    data <- as.data.frame(data)
  } else if (!inherits(data, "mids")) {
    stop("'data' must be either a data.frame or mids object.")
  }

  check_inputs_adjustedcif(data=data, variable=variable, ev_time=ev_time,
                           event=event, cause=cause, method=method,
                           conf_int=conf_int, conf_level=conf_level,
                           times=times, bootstrap=bootstrap,
                           n_boot=n_boot, na.action=na.action,
                           clean_data=clean_data, ...)

  # get required packages
  three_dots <- list(...)
  load_needed_packages(method=method, kind="cif",
                       treatment_model=three_dots$treatment_model,
                       censoring_vars=three_dots$censoring_vars)

  if (inherits(data, "mids")) {

    # levels of the group variable
    if (is.numeric(data$data[, variable])) {
      levs <- unique(data$data[, variable])
    } else {
      levs <- levels(data$data[, variable])
    }

    # transform to long format
    mids <- mice::complete(data, action="long", include=FALSE)

    # get event specific times
    if (is.null(times)) {
      times <- sort(unique(mids[, ev_time][mids[, event]>=1]))

      # add zero if not already in there
      if (!0 %in% times & method!="tmle") {
        times <- c(0, times)
      }
    }

    # get additional arguments
    args <- list(...)

    # extract outcome models
    outcome_models <- args$outcome_model$analyses
    args$outcome_model <- NULL

    # extract treatment models
    if (inherits(args$treatment_model$analyses[[1]],
                 c("glm", "lm", "multinom"))) {
      treatment_models <- args$treatment_model$analyses
      args$treatment_model <- NULL
    } else if (inherits(args$treatment_model, "formula")) {
      treatment_models <- rep(list(args$treatment_model), max(mids$.imp))
      args$treatment_model <- NULL
    } else {
      treatment_models <- NULL
    }

    # extract censoring models
    censoring_models <- args$censoring_model$analyses
    args$censoring_model <- NULL

    # call adjustedcif once for each multiply imputed dataset
    out <- vector(mode="list", length=max(mids$.imp))
    for (i in seq_len(max(mids$.imp))) {

      imp_data <- mids[mids$.imp==i, ]

      # NOTE: need to add the data to the model object or ate() fails
      if (!is.null(treatment_models) &&
          inherits(treatment_models[[i]], "glm")) {
        treatment_models[[i]]$data <- imp_data
      }

      args2 <- list(variable=variable, ev_time=ev_time, event=event,
                    cause=cause, method=method, conf_int=conf_int,
                    conf_level=conf_level, times=times,
                    bootstrap=bootstrap, n_boot=n_boot, n_cores=n_cores,
                    na.action="na.pass", clean_data=clean_data)
      args2 <- c(args2, args)
      args2$data <- imp_data
      args2$outcome_model <- outcome_models[[i]]
      args2$treatment_model <- treatment_models[[i]]
      args2$censoring_model <- censoring_models[[i]]

      out[[i]] <- do.call(adjustedcif, args=args2)

    }

    # pool results
    dats <- vector(mode="list", length=length(out))
    boot_dats <- vector(mode="list", length=length(out))
    for (i in seq_len(length(out))) {

      # direct estimate
      dat <- out[[i]]$adj
      dat$.imp <- i
      dats[[i]] <- dat

      # bootstrap estimate
      boot_dat <- out[[i]]$boot_adj
      boot_dat$.imp <- i
      boot_dats[[i]] <- boot_dat

    }
    dats <- dplyr::bind_rows(dats)
    boot_dats <- dplyr::bind_rows(boot_dats)

    if (conf_int) {
      # use Rubins Rule
      plotdata <- dats %>%
        dplyr::group_by(., time, group) %>%
        dplyr::summarise(cif_est=mean(cif, na.rm=mi_extrapolation),
                         var_w = mean(se^2, na.rm = mi_extrapolation),
                         # Estimated between imputation variance
                         var_b = stats::var(cif, na.rm = mi_extrapolation),
                         # Number of imputed datasets
                         B = dplyr::n(),
                         # Estimated total variance
                         var_t = var_w + var_b + var_b/B,
                         se = sqrt(var_t),
                         .groups="drop_last") %>%
        # dplyr::select(-var_w, -var_b, -B, -var_t) %>%
        dplyr::select(-B) %>%
        dplyr::rename(cif = cif_est)

      plotdata <- as.data.frame(plotdata)

      # re-calculate confidence intervals using pooled se
      surv_ci <- confint_surv(surv=plotdata$cif,
                              se=plotdata$se,
                              conf_level=conf_level,
                              conf_type="plain")
      plotdata$ci_lower <- surv_ci$left
      plotdata$ci_upper <- surv_ci$right
    } else {
      plotdata <- dats %>%
        dplyr::group_by(., time, group) %>%
        dplyr::summarise(cif=mean(cif, na.rm=mi_extrapolation),
                         .groups="drop_last")
      plotdata <- as.data.frame(plotdata)
    }

    # if estimated CIFs beyond the maximum observed survival time
    # exist, remove them if specified
    if (!mi_extrapolation) {
      max_t_group <- max_observed_time(mids=data, variable=variable,
                                       ev_time=ev_time, event=event,
                                       cause=cause, levs=levs,
                                       method=method, type="cif")
      plotdata <- merge(plotdata, max_t_group, by="group", all.x=TRUE)
      plotdata <- plotdata[plotdata$time <= plotdata$max_t,]
      plotdata$max_t <- NULL
      plotdata <- plotdata[order(plotdata$group, plotdata$time), ]
    }

    if (force_bounds) {
      plotdata <- force_bounds_est(plotdata)
    }

    if (iso_reg) {
      plotdata <- iso_reg_est(plotdata)
    }

    # output object
    out_obj <- list(mids_analyses=out,
                    adj=plotdata,
                    data=data$data,
                    method=method,
                    categorical=ifelse(length(levs)>2, TRUE, FALSE),
                    call=match.call())

    if (bootstrap) {

      plotdata_boot <- boot_dats %>%
        dplyr::group_by(., time, group) %>%
        dplyr::summarise(cif=mean(cif),
                         se=mean(se),
                         .groups="drop_last")
      plotdata_boot <- as.data.frame(plotdata_boot)

      # re-calculate confidence intervals using pooled se
      surv_ci <- confint_surv(surv=plotdata_boot$cif,
                              se=plotdata_boot$se,
                              conf_level=conf_level,
                              conf_type="plain")
      plotdata_boot$ci_lower <- surv_ci$left
      plotdata_boot$ci_upper <- surv_ci$right

      out_obj$boot_adj <- plotdata_boot

    }

    class(out_obj) <- "adjustedcif"
    return(out_obj)

  ## normal method using a single data.frame
  } else {

    # only keep needed covariates
    if (clean_data) {
      data <- remove_unnecessary_covars(data=data, variable=variable,
                                        method=method, ev_time=ev_time,
                                        event=event, ...)
    }

    # cant move to input checks because it has to be called after
    # removal of useless covariates
    if (anyNA(data) && is.numeric(three_dots$treatment_model)) {
      stop("Weights cannot be supplied directly to the 'treatment_model'",
           " argument if there are missing values in relevant",
           " columns of 'data'.")
    }

    # perform na.action
    if (is.function(na.action)) {
      data <- na.action(data)
    } else {
      na.action <- get(na.action)
      data <- na.action(data)
    }

    # edge case: there is no data left after removals
    if (nrow(data)==0) {
      stop("There is no non-missing data left after call to 'na.action'.")
    }

    # define those to remove Notes in devtools::check()
    . <- i <- time <- group <- cif_b <- cif <- se <- NULL

    # get event specific times
    if (is.null(times) & !bootstrap & method=="aalen_johansen") {
      times <- NULL
    } else if (is.null(times)) {
      times <- sort(unique(data[, ev_time][data[, event]>=1]))

      # add zero if not already in there
      if (!0 %in% times & method!="tmle") {
        times <- c(0, times)
      }
    }

    # levels of the group variable
    if (is.numeric(data[, variable])) {
      levs <- unique(data[, variable])
    } else {
      levs <- levels(data[, variable])
    }

    # get relevant cif_method function
    cif_fun <- get(paste0("cif_", method))

    # bootstrap the whole procedure, can be useful to get se, p-values
    if (bootstrap) {

      if (n_cores > 1) {
        # needed packages for parallel processing
        requireNamespace("parallel")
        requireNamespace("doRNG")
        requireNamespace("doParallel")
        requireNamespace("foreach")

        # initialize clusters
        cl <- parallel::makeCluster(n_cores, outfile="")
        doParallel::registerDoParallel(cl)
        pkgs <- (.packages())
        export_objs <- c("get_iptw_weights", "read_from_fun",
                         "adjustedcif_boot", "trim_weights",
                         "geese_predictions", "load_needed_packages",
                         "specific_times", "cif_g_comp")

        boot_out <- foreach::foreach(i=1:n_boot, .packages=pkgs,
                                    .export=export_objs) %dorng% {

          adjustedcif_boot(data=data, variable=variable, ev_time=ev_time,
                           event=event, method=method, times=times, i=i,
                           cif_fun=cif_fun, cause=cause,
                           na.action=na.action, ...)
        }
        parallel::stopCluster(cl)

      } else {

        boot_out <- vector(mode="list", length=n_boot)
        for (i in seq_len(n_boot)) {
          boot_out[[i]] <- adjustedcif_boot(data=data, variable=variable,
                                            ev_time=ev_time, event=event,
                                            method=method,
                                            times=times, i=i, cause=cause,
                                            cif_fun=cif_fun,
                                            na.action=na.action, ...)
        }
      }

      # transform into data.frame
      boot_data <- as.data.frame(dplyr::bind_rows(boot_out))

      # keep factor ordering the same
      boot_data$group <- factor(boot_data$group, levels=levs)

      # calculate some statistics
      boot_stats <- boot_data %>%
        dplyr::group_by(., time, group) %>%
        dplyr::summarise(boot_cif=mean(cif, na.rm=TRUE),
                         se=stats::sd(cif, na.rm=TRUE),
                         ci_lower=stats::quantile(cif,
                                                  probs=(1-conf_level)/2,
                                                  na.rm=TRUE),
                         ci_upper=stats::quantile(cif,
                                                  probs=1-((1-conf_level)/2),
                                                  na.rm=TRUE),
                         n_boot=sum(!is.na(cif)),
                         .groups="drop_last")
      boot_stats$group <- factor(boot_stats$group, levels=levs)
    }

    # core of the function
    args <- list(data=data, variable=variable, ev_time=ev_time,
                 event=event, conf_int=conf_int, conf_level=conf_level,
                 times=times, cause=cause, ...)
    method_results <- R.utils::doCall(cif_fun, args=args,
                                      .ignoreUnusedArgs=FALSE)
    plotdata <- method_results$plotdata

    # keep factor ordering the same
    plotdata$group <- factor(plotdata$group, levels=levs)

    if (force_bounds) {
      plotdata <- force_bounds_est(plotdata)
    }

    if (iso_reg) {
      plotdata <- iso_reg_est(plotdata)
    }

    out <- list(adj=plotdata,
                data=data,
                method=method,
                categorical=ifelse(length(levs)>2, TRUE, FALSE),
                call=match.call())

    if (bootstrap) {
      out$boot_data <- boot_data

      # relevant estimates
      plotdata_temp <- dplyr::select(plotdata, c("time", "group", "cif"))
      plotdata_temp <- as.data.frame(plotdata_temp)
      boot_stats <- as.data.frame(boot_stats)

      # order both data.frames
      plotdata_temp <- plotdata_temp[order(plotdata_temp$group,
                                           plotdata_temp$time),]
      boot_stats <- boot_stats[order(boot_stats$group,
                                     boot_stats$time),]

      # put together
      boot_stats$cif <- plotdata_temp$cif
      out$boot_adj <- boot_stats
    }

    # add method-specific objects to output
    method_results$plotdata <- NULL
    out <- c(out, method_results)

    class(out) <- "adjustedcif"
    return(out)
  }
}

## perform one bootstrap iteration
adjustedcif_boot <- function(data, variable, ev_time, event, cause, method,
                             times, i, cif_fun, na.action, ...) {

  # draw sample
  indices <- sample(x=rownames(data), size=nrow(data), replace=TRUE)
  boot_samp <- data[indices, ]

  # perform na.action
  boot_samp <- na.action(boot_samp)

  # IMPORTANT: keeps SL in tmle methods from failing
  row.names(boot_samp) <- seq_len(nrow(data))

  # update models/recalculate weights using bootstrap sample
  pass_args <- list(...)
  if (method %in% c("direct", "aiptw") &
      !inherits(pass_args$outcome_model, "formula")) {
    pass_args$outcome_model <- stats::update(pass_args$outcome_model,
                                             data=boot_samp)
  }

  if ((method %in% c("iptw", "iptw_pseudo", "aiptw", "aiptw_pseudo")) &&
      inherits(pass_args$treatment_model, c("glm", "multinom"))) {
    pass_args$treatment_model <- stats::update(pass_args$treatment_model,
                                               data=boot_samp, trace=FALSE)
  }

  if ((method %in% c("iptw", "direct", "aiptw")) &&
      inherits(pass_args$censoring_model, "coxph")) {
    pass_args$censoring_model <- stats::update(pass_args$censoring_model,
                                               data=boot_samp)
  }

  # call cif_method with correct arguments
  args <- list(data=boot_samp, variable=variable, ev_time=ev_time,
               event=event, conf_int=FALSE, conf_level=0.95, times=times,
               cause=cause)
  args <- c(args, pass_args)

  method_results <- R.utils::doCall(cif_fun, args=args, .ignoreUnusedArgs=FALSE)
  adj_boot <- method_results$plotdata
  adj_boot$boot <- i

  return(adj_boot)

}

## S3 summary method for adjustedcif objects
#' @export
summary.adjustedcif <- function(object, ...) {

  if (object$method=="direct") {
    method_name <- "Direct Standardization"
  } else if (object$method=="direct_pseudo") {
    method_name <- "Direct Standardization: Pseudo-Values"
  } else if (object$method=="iptw") {
    method_name <- "Inverse Probability of Treatment Weighting"
  } else if (object$method=="iptw_pseudo") {
    method_name <- "Inverse Probability of Treatment Weighting: Pseudo-Values"
  } else if (object$method=="matching") {
    method_name <- "Propensity Score Matching"
  } else if (object$method=="aiptw") {
    method_name <- "Augmented Inverse Probability of Treatment Weighting"
  } else if (object$method=="aiptw_pseudo") {
    method_name <- paste0("Augmented Inverse Probability of Treatment",
                          " Weighting: Pseudo-Values")
  } else if (object$method=="aalen_johansen") {
    method_name <- "Aalen-Johansen Estimator"
  } else if (object$method=="tmle") {
    method_name <- "Targeted Maximum Likelihood Estimator"
  }

  times_str <- ifelse(is.null(object$call$times), "Event-Specific Times",
                      "User-Supplied Points in Time")

  if (object$method=="aalen_johansen") {
    cat("Unadjusted Cumulative Incidences \n")
  } else {
    cat("Confounder Adjusted Cumulative Incidences \n")
  }
  cat("   - Cause of Interest: ", object$call$cause, "\n", sep="")
  cat("   - Method: ", method_name, "\n", sep="")
  cat("   - Times: ", times_str, "\n", sep="")

  if (is.null(object$call$bootstrap) ||
      !as.logical(as.character(object$call$bootstrap))) {
    cat("   - Bootstrapping: Not Done\n", sep="")
  } else {
    n_boot <- ifelse(is.null(object$call$n_boot), 500, object$call$n_boot)
    cat("   - Bootstrapping: Performed with ", n_boot,
        " Replications\n", sep="")
  }

  if (is.null(object$call$conf_int) ||
      !as.logical(as.character(object$call$conf_int))) {
    cat("   - Approximate CI: Not Calculated\n", sep="")
  } else {
    conf_level <- ifelse(is.null(object$call$conf_level), 0.95,
                         object$call$conf_level)
    cat("   - Approximate CI: Calculated with a Confidence level of ",
        conf_level, "\n", sep="")
  }

  if (is.null(object$mids_analyses)) {
    cat("   - Using a single dataset\n")
  } else {
    cat("   - Using multiply imputed dataset\n")
  }
}

## S3 print method for adjustedcif objects
#' @export
print.adjustedcif <- function(x, ...) {
  summary(x, ...)
}

## S3 print method for adjustedcif.method objects
#' @export
print.adjustedcif.method <- function(x, ...) {
  print(x$plotdata, ...)
}

## S3 summary method for adjustedcif.method objects
#' @export
summary.adjustedcif.method <- function(object, ...) {
  summary(object$plotdata, ...)
}
