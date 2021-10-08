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

## Main function of the package. Is basically a wrapper around
## all other functions, offering additional high level stuff
#' @importFrom dplyr %>%
#' @importFrom doRNG %dorng%
#' @importFrom foreach %dopar%
#' @export
adjustedsurv <- function(data, variable, ev_time, event, method, conf_int=FALSE,
                         conf_level=0.95, times=NULL, bootstrap=FALSE,
                         n_boot=500, n_cores=1,
                         na.action=options("na.action")[[1]],
                         ...) {

  check_inputs_adjustedsurv(data=data, variable=variable,
                            ev_time=ev_time, event=event, method=method,
                            conf_int=conf_int, conf_level=conf_level,
                            times=times, bootstrap=bootstrap,
                            n_boot=n_boot, na.action=na.action, ...)

  # get required packages
  three_dots <- list(...)
  load_needed_packages(method=method, kind="surv",
                       treatment_model=three_dots$treatment_model,
                       censoring_vars=three_dots$censoring_vars)
  rm(three_dots)

  ## using multiple imputation
  if (inherits(data, "mids")) {

    # get event specific times
    times_input <- times
    if (is.null(times)) {
      times <- sort(unique(data$data[, ev_time][data$data[, event]==1]))

      # add zero if not already in there
      if (!0 %in% times) {
        times <- c(0, times)
      }
    }

    # levels of the group variable
    if (is.numeric(data$data[,variable])) {
      levs <- unique(data$data[,variable])
    } else {
      levs <- levels(data$data[,variable])
    }

    # transform to long format
    mids <- mice::complete(data, action="long", include=FALSE)

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
    } else if (is.numeric(args$treatment_model)) {
      stop("Supplying weights or propensity scores directly is not allowed",
           " when using multiple imputation.")
    } else {
      treatment_models <- NULL
    }

    # extract censoring models
    censoring_models <- args$censoring_model$analyses
    args$censoring_model <- NULL

    # call adjustedsurv once for each multiply imputed dataset
    out <- vector(mode="list", length=max(mids$.imp))
    for (i in seq_len(max(mids$.imp))) {

      imp_data <- mids[mids$.imp==i,]

      # NOTE: need to add the data to the model object or ate() fails
      if (!is.null(treatment_models) &
          inherits(treatment_models[[i]], "glm")) {
        treatment_models[[i]]$data <- imp_data
      }

      args2 <- c(variable=variable, ev_time=ev_time,
                 event=event, method=method, conf_int=conf_int,
                 conf_level=conf_level, times=times,
                 bootstrap=bootstrap, n_boot=n_boot, n_cores=n_cores,
                 na.action="na.pass", args)
      args2$data <- imp_data
      args2$outcome_model <- outcome_models[[i]]
      args2$treatment_model <- treatment_models[[i]]
      args2$censoring_model <- censoring_models[[i]]

      out[[i]] <- do.call(adjustedsurv, args=args2)

    }

    # pool results
    dats <- vector(mode="list", length=length(out))
    boot_dats <- vector(mode="list", length=length(out))
    for (i in seq_len(length(out))) {

      # direct estimate
      dat <- out[[i]]$adjsurv
      dat$.imp <- i
      dats[[i]] <- dat

      # bootstrap estimate
      boot_dat <- out[[i]]$boot_adjsurv
      boot_dat$.imp <- i
      boot_dats[[i]] <- boot_dat

    }
    dats <- dplyr::bind_rows(dats)
    boot_dats <- dplyr::bind_rows(boot_dats)

    if (conf_int) {
      # use Rubins Rule
      plotdata <- dats %>%
        dplyr::group_by(., time, group) %>%
        dplyr::summarise(surv=mean(surv),
                         se=mean(se),
                         .groups="drop_last")
      plotdata <- as.data.frame(plotdata)

      # re-calculate confidence intervals using pooled se
      surv_ci <- confint_surv(surv=plotdata$surv,
                              se=plotdata$se,
                              conf_level=conf_level,
                              conf_type="plain")
      plotdata$ci_lower <- surv_ci$left
      plotdata$ci_upper <- surv_ci$right
    } else {
      plotdata <- dats %>%
        dplyr::group_by(., time, group) %>%
        dplyr::summarise(surv=mean(surv),
                         .groups="drop_last")
      plotdata <- as.data.frame(plotdata)
    }

    # output object
    out_obj <- list(mids_analyses=out,
                    adjsurv=plotdata,
                    data=data$data,
                    method=method,
                    categorical=ifelse(length(levs)>2, TRUE, FALSE),
                    call=match.call())

    if (bootstrap) {

      plotdata_boot <- boot_dats %>%
        dplyr::group_by(., time, group) %>%
        dplyr::summarise(surv=mean(surv),
                         se=mean(se),
                         .groups="drop_last")
      plotdata_boot <- as.data.frame(plotdata_boot)

      # re-calculate confidence intervals using pooled se
      surv_ci <- confint_surv(surv=plotdata_boot$surv,
                              se=plotdata_boot$se,
                              conf_level=conf_level,
                              conf_type="plain")
      plotdata_boot$ci_lower <- surv_ci$left
      plotdata_boot$ci_upper <- surv_ci$right

      out_obj$boot_adjsurv <- plotdata_boot

    }

  class(out_obj) <- "adjustedsurv"
  return(out_obj)

  ## normal method using a single data.frame
  } else {

    # only keep needed covariates
    data <- remove_unnecessary_covars(data=data, variable=variable,
                                      method=method, ev_time=ev_time,
                                      event=event, ...)

    # perform na.action
    if (is.function(na.action)) {
      data <- na.action(data)
    } else {
      na.action <- get(na.action)
      data <- na.action(data)
    }

    # define those to remove Notes in devtools::check()
    . <- i <- time <- group <- surv_b <- surv <- se <- NULL

    # get event specific times
    times_input <- times
    if (is.null(times) & method %in% c("km", "iptw_km")) {
      times <- NULL
    } else if (is.null(times)) {
      times <- sort(unique(data[, ev_time][data[, event]==1]))

      # add zero if not already in there
      if (!0 %in% times) {
        times <- c(0, times)
      }
    }

    # levels of the group variable
    if (is.numeric(data[,variable])) {
      levs <- unique(data[,variable])
    } else {
      levs <- levels(data[,variable])
    }

    # get relevant surv_method function
    surv_fun <- get(paste0("surv_", method))

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
        pkgs <- c("adjustedCurves", "survival")
        export_objs <- c("get_iptw_weights", "read_from_step_function",
                         "multi_result_class", "adjustedsurv_boot",
                         "trim_weights", "calc_pseudo_surv",
                         "geese_predictions", "load_needed_packages",
                         "specific_times")

        boot_out <- foreach::foreach(i=1:n_boot, .packages=pkgs,
                                     .export=export_objs) %dorng% {

          adjustedsurv_boot(data=data, variable=variable, ev_time=ev_time,
                            event=event, method=method, times_input=times_input,
                            times=times, i=i, surv_fun=surv_fun,
                            levs=levs, na.action=na.action, ...)
                                     }
        parallel::stopCluster(cl)

      } else {

        boot_out <- vector(mode="list", length=n_boot)
        for (i in seq_len(n_boot)) {
          boot_out[[i]] <- adjustedsurv_boot(data=data, variable=variable,
                                             ev_time=ev_time, event=event,
                                             method=method,
                                             times_input=times_input,
                                             times=times, i=i,
                                             surv_fun=surv_fun, levs=levs,
                                             na.action=na.action, ...)
        }
      }

      # transform into data.frame
      boot_data <- as.data.frame(dplyr::bind_rows(boot_out))

      # keep factor ordering the same
      boot_data$group <- factor(boot_data$group, levels=levs)
      colnames(boot_data) <- c("time", "surv_b", "group", "boot")

      # calculate some statistics
      boot_stats <- boot_data %>%
        dplyr::group_by(., time, group) %>%
        dplyr::summarise(surv=mean(surv_b, na.rm=TRUE),
                         se=stats::sd(surv_b, na.rm=TRUE),
                         ci_lower=stats::quantile(surv_b,
                                                  probs=(1-conf_level)/2,
                                                  na.rm=TRUE),
                         ci_upper=stats::quantile(surv_b,
                                                  probs=1-((1-conf_level)/2),
                                                  na.rm=TRUE),
                         n_boot=sum(!is.na(surv_b)),
                         .groups="drop_last")
      boot_stats$group <- factor(boot_stats$group, levels=levs)

      # get old names back
      colnames(boot_data) <- c("time", "surv", "group", "boot")
    }

    # core of the function
    args <- list(data=data, variable=variable, ev_time=ev_time,
                 event=event, conf_int=conf_int, conf_level=conf_level,
                 times=times, ...)
    method_results <- R.utils::doCall(surv_fun, args=args)
    plotdata <- method_results$plotdata

    # keep factor levels in same order as data
    plotdata$group <- factor(plotdata$group, levels=levs)

    out <- list(adjsurv=plotdata,
                data=data,
                method=method,
                categorical=ifelse(length(levs) > 2, TRUE, FALSE),
                call=match.call())

    if (bootstrap) {
      out$boot_data <- boot_data
      out$boot_adjsurv <- as.data.frame(boot_stats)
    }

    # add method-specific objects to output
    method_results$plotdata <- NULL
    out <- c(out, method_results)

    class(out) <- "adjustedsurv"
    return(out)

  }
}

## perform one bootstrap iteration
adjustedsurv_boot <- function(data, variable, ev_time, event, method,
                              times_input, times, i, surv_fun, levs,
                              na.action, ...) {

  # get required packages
  three_dots <- list(...)
  load_needed_packages(method=method, kind="surv",
                       treatment_model=three_dots$treatment_model,
                       censoring_vars=three_dots$censoring_vars)
  rm(three_dots)

  # draw sample
  indices <- sample(x=rownames(data), size=nrow(data), replace=TRUE)
  boot_samp <- data[indices,]

  # perform na.action
  if (is.function(na.action)) {
    boot_samp <- na.action(boot_samp)
  } else {
    na.action <- get(na.action)
    boot_samp <- na.action(boot_samp)
  }

  # IMPORTANT: keeps SL in tmle methods from failing
  row.names(boot_samp) <- seq_len(nrow(data))

  # update models/recalculate weights using bootstrap sample
  pass_args <- list(...)
  if ((method %in% c("direct", "aiptw")) &
      !inherits(pass_args$outcome_model, "formula")) {
    pass_args$outcome_model <- stats::update(pass_args$outcome_model,
                                             data=boot_samp)
  }

  if ((method %in% c("iptw_km", "iptw_cox", "iptw_pseudo", "aiptw",
                    "aiptw_pseudo")) &&
      (inherits(pass_args$treatment_model, c("glm", "multinom")))) {
    pass_args$treatment_model <- stats::update(pass_args$treatment_model,
                                               data=boot_samp, trace=FALSE)
  }

  if ((method %in% c("direct", "aiptw")) &&
      inherits(pass_args$censoring_model, "coxph")) {
    pass_args$censoring_model <- stats::update(pass_args$censoring_model,
                                               data=boot_samp)
  }

  # call surv_method with correct arguments
  args <- list(data=boot_samp, variable=variable, ev_time=ev_time,
               event=event, conf_int=FALSE, conf_level=0.95, times=times)
  args <- c(args, pass_args)

  method_results <- R.utils::doCall(surv_fun, args=args)
  adjsurv_boot <- method_results$plotdata
  adjsurv_boot$boot <- i

  return(adjsurv_boot)
}

## plot the survival curves
#' @importFrom rlang .data
#' @export
plot.adjustedsurv <- function(x, draw_ci=FALSE, max_t=Inf,
                              iso_reg=FALSE, force_bounds=FALSE,
                              use_boot=FALSE, color=TRUE,
                              linetype=FALSE, facet=FALSE,
                              line_size=1, xlab="Time",
                              ylab="Adjusted Survival Probability",
                              title=NULL, legend.title="Group",
                              legend.position="right",
                              gg_theme=ggplot2::theme_classic(),
                              ylim=NULL, custom_colors=NULL,
                              custom_linetypes=NULL,
                              ci_draw_alpha=0.4, steps=TRUE,
                              median_surv_lines=FALSE, median_surv_size=0.5,
                              median_surv_linetype="dashed",
                              median_surv_color="black",
                              censoring_ind=FALSE, censoring_ind_width=NULL,
                              censoring_ind_size=0.5, ...) {
  requireNamespace("ggplot2")

  if (use_boot & is.null(x$boot_adjsurv)) {
    warning("Cannot use bootstrapped estimates as they were not estimated.",
            " Need bootstrap=TRUE in adjustedsurv() call.")
    draw_ci <- FALSE
    plotdata <- x$adjsurv
  } else if (use_boot) {
    plotdata <- x$boot_adjsurv
  } else {
    plotdata <- x$adjsurv
  }
  plotdata$group <- factor(plotdata$group)

  # shortcut to only show curves up to a certain time
  plotdata <- plotdata[which(plotdata$time <= max_t),]

  # in some methods estimates can be outside the 0, 1 bounds,
  # if specified set those to 0 or 1 respectively
  if (force_bounds) {
    plotdata <- within(plotdata, {
      surv <- ifelse(surv < 0, 0, surv)
      surv <- ifelse(surv > 1, 1, surv)
    })
  }

  # apply isotonic regression if specified
  if (iso_reg) {
    for (lev in levels(plotdata$group)) {
      temp <- plotdata[plotdata$group==lev,]
      # to surv estimates
      new <- rev(stats::isoreg(rev(temp$surv))$yf)
      plotdata$surv[plotdata$group==lev] <- new
      # shift confidence intervals accordingly
      if (draw_ci & "ci_lower" %in% colnames(temp)) {
        diff <- temp$surv - new

        plotdata$ci_lower[plotdata$group==lev] <- temp$ci_lower - diff
        plotdata$ci_upper[plotdata$group==lev] <- temp$ci_upper - diff
      }
    }
  }

  ## The main plot
  mapping <- ggplot2::aes(x=.data$time, y=.data$surv, color=.data$group,
                          linetype=.data$group, group=.data$group)

  if (!linetype) {
    mapping$linetype <- NULL
  }
  if (!color) {
    mapping$colour <- NULL
  }

  p <- ggplot2::ggplot(plotdata, mapping)

  if (steps) {
    p <- p + ggplot2::geom_step(size=line_size)
  } else {
    p <- p + ggplot2::geom_line(size=line_size)
  }

  p <- p + gg_theme +
    ggplot2::labs(x=xlab, y=ylab, color=legend.title,
                  linetype=legend.title, fill=legend.title) +
    ggplot2::theme(legend.position=legend.position)

  if (facet) {
    p <- p + ggplot2::facet_wrap(~group)
  }
  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }
  if (!is.null(ylim)) {
    p <- p + ggplot2::ylim(ylim)
  }
  if (!is.null(custom_colors)) {
    p <- p + ggplot2::scale_colour_manual(values=custom_colors)
  }
  if (!is.null(custom_linetypes)) {
    p <- p + ggplot2::scale_linetype_manual(values=custom_linetypes)
  }
  if (!is.null(custom_colors)) {
    p <- p + ggplot2::scale_fill_manual(values=custom_colors)
  }

  ## Censoring indicators
  if (censoring_ind) {

    if (is.null(censoring_ind_width)) {

      if (is.null(ylim)) {
        ystart <- 1 - ggplot2::layer_scales(p)$y$range$range[1]
      } else {
        ystart <- 1 - ylim[1]
      }

      censoring_ind_width <- ystart * 0.05

    }

    if (is.factor(plotdata$group)) {
      levs <- levels(plotdata$group)
    } else {
      levs <- unique(plotdata$group)
    }

    # calculate needed data
    cens_dat <- vector(mode="list", length=length(levs))
    for (i in seq_len(length(levs))) {

      x$data <- x$data[which(x$data$time <= max_t),]
      cens_times <- sort(unique(x$data[, x$call$ev_time][
        x$data[, x$call$event]==0 & x$data[, x$call$variable]==levs[i]]))
      adjsurv_temp <- plotdata[plotdata$group==levs[i], ]
      cens_surv <- vapply(cens_times, read_from_step_function,
                          step_data=adjsurv_temp, FUN.VALUE=numeric(1))
      cens_dat[[i]] <- data.frame(time=cens_times, surv=cens_surv,
                                  group=levs[i])

    }
    cens_dat <- dplyr::bind_rows(cens_dat)
    cens_dat <- cens_dat[!is.na(cens_dat$surv) ,]

    cens_map <- ggplot2::aes(x=.data$time,
                             y=.data$surv-(censoring_ind_width/2),
                             xend=.data$time,
                             yend=.data$surv+(censoring_ind_width/2),
                             group=.data$group,
                             color=.data$group,
                             linetype=.data$group)
    if (!color) {
      cens_map$colour <- NULL
    }
    if (!linetype) {
      cens_map$linetype <- NULL
    }

    p <- p + ggplot2::geom_segment(data=cens_dat, cens_map,
                                   size=censoring_ind_size)

  }

  ## Confidence intervals
  if (draw_ci & !"ci_lower" %in% colnames(plotdata)) {
    warning("Cannot draw confidence intervals. Either set 'conf_int=TRUE' in",
            " 'adjustedsurv()' call or use bootstrap estimates.")
  }

  if ((draw_ci & "ci_lower" %in% colnames(plotdata)) & steps) {
    ci_map <- ggplot2::aes(ymin=.data$ci_lower,
                           ymax=.data$ci_upper,
                           group=.data$group,
                           fill=.data$group,
                           x=.data$time,
                           y=.data$surv)

    if (!color) {
      ci_map$fill <- NULL
    }

    p <- p + pammtools::geom_stepribbon(ci_map, alpha=ci_draw_alpha,
                                        inherit.aes=FALSE)
  } else if (draw_ci & "ci_lower" %in% colnames(plotdata)) {
    ci_map <- ggplot2::aes(ymin=.data$ci_lower,
                           ymax=.data$ci_upper,
                           group=.data$group,
                           fill=.data$group,
                           x=.data$time,
                           y=.data$surv)

    if (!color) {
      ci_map$fill <- NULL
    }

    p <- p + ggplot2::geom_ribbon(ci_map, alpha=ci_draw_alpha,
                                  inherit.aes=FALSE)
  }

  ## Median Survival indicators
  if (median_surv_lines) {

    # calculate median survival and add other needed values
    fake_adjsurv <- x
    fake_adjsurv$adjsurv <- plotdata
    median_surv <- adjusted_median_survival(fake_adjsurv, use_boot=FALSE,
                                            verbose=FALSE)
    median_surv$y <- 0.5
    # set to NA if not in plot
    median_surv$median_surv[median_surv$median_surv > max_t] <- NA

    if (is.null(ylim)) {
      median_surv$yend <- ggplot2::layer_scales(p)$y$range$range[1]
    } else {
      median_surv$yend <- ylim[1]
    }

    # remove if missing
    median_surv <- median_surv[!is.na(median_surv$median_surv),]

    if (sum(is.na(median_surv$median_surv)) < nrow(median_surv)) {

      # draw line on surv_p = 0.5 until it hits the last curve
      p <- p + ggplot2::geom_segment(ggplot2::aes(x=0,
                                              xend=max(median_surv$median_surv),
                                              y=0.5,
                                              yend=0.5),
                                     inherit.aes=FALSE,
                                     linetype=median_surv_linetype,
                                     size=median_surv_size,
                                     color=median_surv_color)
      # draw indicator lines from middle to bottom
      p <- p + ggplot2::geom_segment(ggplot2::aes(x=.data$median_surv,
                                                  xend=.data$median_surv,
                                                  y=0.5,
                                                  yend=.data$yend),
                                     inherit.aes=FALSE,
                                     linetype=median_surv_linetype,
                                     size=median_surv_size,
                                     color=median_surv_color,
                                     data=median_surv)
    }
  }
  return(p)
}

## S3 print method for adjustedsurv objects
#' @export
print.adjustedsurv <- function(x, ...) {

  if (x$method=="direct") {
    method_name <- "Direct Standardization"
  } else if (x$method=="direct_pseudo") {
    method_name <- "Direct Standardization: Pseudo-Values"
  } else if (x$method=="iptw_km") {
    method_name <- "Inverse Probability of Treatment Weighting: Kaplan-Meier"
  } else if (x$method=="iptw_cox") {
    method_name <- "Inverse Probability of Treatment Weighting: Cox-Regression"
  } else if (x$method=="iptw_pseudo") {
    method_name <- "Inverse Probability of Treatment Weighting: Pseudo-Values"
  } else if (x$method=="matching") {
    method_name <- "Propensity Score Matching"
  } else if (x$method=="emp_lik") {
    method_name <- "Empirical Likelihood Estimation"
  } else if (x$method=="aiptw") {
    method_name <- "Augmented Inverse Probability of Treatment Weighting"
  } else if (x$method=="aiptw_pseudo") {
    method_name <- paste0("Augmented Inverse Probability of Treatment",
                          " Weighting: Pseudo-Values")
  } else if (x$method=="tmle") {
    method_name <- "Targeted Maximum Likelihood Estimation"
  } else if (x$method=="ostmle") {
    method_name <- "One-Step Targeted Maximum Likelihood Estimation"
  } else if (x$method=="km") {
    method_name <- "Kaplan-Meier Estimator"
  }

  times_str <- ifelse(is.null(x$call$times), "Event-Specific Times",
                      "User-Supplied Points in Time")

  if (x$method=="km") {
    cat("Unadjusted Survival Probabilities \n")
  } else {
    cat("Confounder Adjusted Survival Probabilities \n")
  }
  cat("   - Method: ", method_name, "\n", sep="")
  cat("   - Times: ", times_str, "\n", sep="")

  if (as.logical(toString(x$call$bootstrap))) {
    cat("   - Bootstrapping: Performed with ", x$call$n_boot,
        " Replications\n", sep="")
  } else {
    cat("   - Bootstrapping: Not none\n", sep="")
  }

  if (is.null(x$call$conf_int) || !as.logical(as.character(x$call$conf_int))) {
    cat("   - Approximate CI: Not calculated\n", sep="")
  } else {
    conf_level <- ifelse(is.null(x$call$conf_level), 0.95, x$call$conf_level)
    cat("   - Approximate CI: Calculated with a confidence level of ",
        conf_level, "\n", sep="")
  }

  if (is.null(x$mids_analyses)) {
    cat("   - Using a single dataset")
  } else {
    cat("   - Using multiply imputed dataset")
  }

}
