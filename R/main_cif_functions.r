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
                        conf_int=F, conf_level=0.95, times=NULL, bootstrap=F,
                        n_boot=500, n_cores=1, na.action=options("na.action")[[1]],
                        ...) {

  check_inputs_adjustedcif(data=data, variable=variable, ev_time=ev_time,
                           event=event, cause=cause, method=method,
                           conf_int=conf_int, conf_level=conf_level,
                           times=times, bootstrap=bootstrap,
                           n_boot=n_boot, na.action=na.action, ...)

  # perform na.action
  if (is.function(na.action)) {
    data <- na.action(data)
  } else {
    na.action <- get(na.action)
    data <- na.action(data)
  }

  # define those to remove Notes in devtools::check()
  . <- i <- time <- group <- cif_b <- NULL

  # get event specific times
  times_input <- times
  if (is.null(times) & method=="aalen_johansen") {
    times <- NULL
  } else if (is.null(times)) {
    times <- sort(unique(data[, ev_time][data[, event]>=1]))

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

  # get relevant cif_method function
  cif_fun <- get(paste0("cif_", method))

  # bootstrap the whole procedure, can be useful to get sd, p-values
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
                       "multi_result_class", "adjustedcif_boot")

      boot_out <- foreach::foreach(i=1:n_boot, .packages=pkgs,
                                   .export=export_objs) %dorng% {

        adjustedcif_boot(data=data, variable=variable, ev_time=ev_time,
                         event=event, method=method, times_input=times_input,
                         times=times, i=i, cif_fun=cif_fun,
                         levs=levs, cause=cause, na.action=na.action, ...)
      }
      parallel::stopCluster(cl)

    } else {

      boot_out <- vector(mode="list", length=n_boot)
      for (i in 1:n_boot) {
        boot_out[[i]] <- adjustedcif_boot(data=data, variable=variable,
                                          ev_time=ev_time, event=event,
                                          method=method, times_input=times_input,
                                          times=times, i=i, cause=cause,
                                          cif_fun=cif_fun, levs=levs,
                                          na.action=na.action, ...)
      }
    }

    # transform into data.frames
    boot_data <- lapply(boot_out, function(x) x$boot_data)
    boot_data_same_t <- lapply(boot_out, function(x) x$boot_data_same_t)

    boot_data <- as.data.frame(dplyr::bind_rows(boot_data))
    boot_data_same_t <- as.data.frame(dplyr::bind_rows(boot_data_same_t))

    # keep factor ordering the same
    boot_data$group <- factor(boot_data$group, levels=levs)
    boot_data_same_t$group <- factor(boot_data_same_t$group, levels=levs)

    # calculate some statistics
    boot_stats <- boot_data_same_t %>%
      dplyr::group_by(., time, group) %>%
      dplyr::summarise(cif=mean(cif_b, na.rm=T),
                       sd=stats::sd(cif_b, na.rm=T),
                       ci_lower=stats::quantile(cif_b,
                                                probs=(1-conf_level)/2,
                                                na.rm=T),
                       ci_upper=stats::quantile(cif_b,
                                                probs=1-((1-conf_level)/2),
                                                na.rm=T),
                       n_boot=sum(!is.na(cif_b)),
                       .groups="drop_last")
    boot_stats$group <- factor(boot_stats$group, levels=levs)
  }

  # core of the function
  args <- list(data=data, variable=variable, ev_time=ev_time,
               event=event, conf_int=conf_int, conf_level=conf_level,
               times=times, cause=cause, ...)
  plotdata <- R.utils::doCall(cif_fun, args=args)

  # keep factor ordering the same
  plotdata$group <- factor(plotdata$group, levels=levs)

  out <- list(adjcif=plotdata,
              method=method,
              categorical=ifelse(length(levs)>2, T, F),
              call=match.call())

  if (bootstrap) {
    out$boot_data <- boot_data
    out$boot_data_same_t <- boot_data_same_t
    out$boot_adjcif <- as.data.frame(boot_stats)
  }

  class(out) <- "adjustedcif"
  return(out)
}

## perform one bootstrap iteration
adjustedcif_boot <- function(data, variable, ev_time, event, cause, method,
                             times_input, times, i, cif_fun, levs,
                             na.action, ...) {

  indices <- sample(x=rownames(data), size=nrow(data), replace=T)
  boot_samp <- data[indices,]

  # perform na.action
  if (is.function(na.action)) {
    boot_samp <- na.action(boot_samp)
  } else {
    na.action <- get(na.action)
    boot_samp <- na.action(boot_samp)
  }

  # IMPORTANT: keeps SL in tmle methods from failing
  row.names(boot_samp) <- 1:nrow(data)

  # if event specific times are used, use event specific times
  # in bootstrapping as well
  if (is.null(times_input) & method=="aalen_johansen") {

  } else if (is.null(times_input)) {
    times_boot <- sort(unique(boot_samp[, ev_time][boot_samp[, event]>=1]))

    if (!0 %in% times_boot) {
      times_boot <- c(0, times_boot)
    }
  } else {
    times_boot <- times
  }

  # update models/recalculate weights using bootstrap sample
  pass_args <- list(...)
  if (method %in% c("direct", "aiptw") &
      !inherits(pass_args$outcome_model, "formula")) {
    pass_args$outcome_model <- stats::update(pass_args$outcome_model,
                                             data=boot_samp)
  }

  if (method %in% c("iptw", "iptw_pseudo", "aiptw", "aiptw_pseudo")) {
    if (inherits(pass_args$treatment_model, "glm") |
        inherits(pass_args$treatment_model, "multinom")) {
      pass_args$treatment_model <- stats::update(pass_args$treatment_model,
                                                 data=boot_samp, trace=F)
    }
  }

  # call cif_method with correct arguments
  args <- list(data=boot_samp, variable=variable, ev_time=ev_time,
               event=event, conf_int=F, conf_level=0.95, times=times,
               cause=cause)
  args <- c(args, pass_args)

  adjcif_boot <- R.utils::doCall(cif_fun, args=args)
  adjcif_boot$boot <- i

  # read from resulting step function at all t in times
  boot_cif <- vector(mode="list", length=length(levs))
  for (j in 1:length(levs)) {

    if (method=="aalen_johansen" & is.null(times)) {
      times <- unique(data[,ev_time][data[,variable]==levs[j]])
    }

    cif_boot <- sapply(times, read_from_step_function,
                       step_data=adjcif_boot[adjcif_boot$group==levs[j],],
                       est="cif")

    dat_temp <- data.frame(time=times,
                           cif_b=cif_boot,
                           group=levs[j],
                           boot=i)
    boot_cif[[j]] <- dat_temp
  }
  boot_cif <- as.data.frame(dplyr::bind_rows(boot_cif))

  # output
  result <- multi_result_class()

  result$boot_data <- adjcif_boot
  result$boot_data_same_t <- boot_cif

  return(result)

}


## plot the cumulative incidence functions
#' @importFrom rlang .data
#' @export
plot.adjustedcif <- function(x, draw_ci=F, max_t=Inf,
                             iso_reg=F, force_bounds=F, use_boot=F,
                             color=T, linetype=F, facet=F,
                             line_size=1, xlab="Time",
                             ylab="Adjusted Cumulative Incidence",
                             title=NULL, legend.title="Group",
                             legend.position="right",
                             ylim=NULL, custom_colors=NULL,
                             custom_linetypes=NULL,
                             ci_draw_alpha=0.4, ...) {

  if (!color & !linetype & !facet) {
    stop("Groups must be distinguished with at least one of 'color',",
         "'linetype' or 'facet'. Can't all be FALSE.")
  }

  if (use_boot & is.null(x$boot_adjcif)) {
    warning("Cannot use bootstrapped estimates as they were not estimated.",
            " Need bootstrap=TRUE in adjustedcif() call.")
    draw_ci <- F
    plotdata <- x$adjcif
  } else if (use_boot) {
    plotdata <- x$boot_adjcif
  } else {
    plotdata <- x$adjcif
  }
  plotdata$group <- factor(plotdata$group)

  # shortcut to only show curves up to a certain time
  plotdata <- plotdata[which(plotdata$time <= max_t),]

  # in some methods estimates can be outside the 0, 1 bounds,
  # if specified set those to 0 or 1 respectively
  if (force_bounds) {
    plotdata <- within(plotdata, {
      cif <- ifelse(cif < 0, 0, cif);
      cif <- ifelse(cif > 1, 1, cif);
    })
  }

  # apply isotonic regression if specified
  if (iso_reg) {
    for (lev in levels(plotdata$group)) {
      cif <- plotdata$cif[plotdata$group==lev]
      new <- stats::isoreg(cif)$yf
      plotdata$cif[plotdata$group==lev] <- new
    }
  }

  mapping <- ggplot2::aes(x=.data$time, y=.data$cif, color=.data$group,
                          linetype=.data$group)

  if (!linetype) {
    mapping$linetype <- NULL
  }
  if (!color) {
    mapping$colour <- NULL
  }

  p <- ggplot2::ggplot(plotdata, mapping) +
    ggplot2::geom_step(size=line_size) +
    ggplot2::theme_bw() +
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

  if (draw_ci & !"ci_lower" %in% colnames(plotdata)) {
    warning("Cannot draw confidence intervals. Either set 'conf_int=TRUE' in",
            " 'adjustedcif()' call or use bootstrap estimates.")
  }

  if (draw_ci & "ci_lower" %in% colnames(plotdata)) {
    p <- p + pammtools::geom_stepribbon(ggplot2::aes(ymin=.data$ci_lower,
                                                     ymax=.data$ci_upper,
                                                     fill=.data$group,
                                                     x=.data$time,
                                                     y=.data$cif),
                                        alpha=ci_draw_alpha, inherit.aes=F)
  }
  return(p)
}

## S3 print method for adjustedcif objects
#' @export
print.adjustedcif <- function(x, ...) {

  if (x$method=="direct") {
    method_name <- "Direct Standardization"
  } else if (x$method=="direct_pseudo") {
    method_name <- "Direct Standardization: Pseudo-Values"
  } else if (x$method=="iptw") {
    method_name <- "Inverse Probability of Treatment Weighting"
  } else if (x$method=="iptw_pseudo") {
    method_name <- "Inverse Probability of Treatment Weighting: Pseudo-Values"
  } else if (x$method=="matching") {
    method_name <- "Propensity Score Matching"
  } else if (x$method=="aiptw") {
    method_name <- "Augmented Inverse Probability of Treatment Weighting"
  } else if (x$method=="aiptw_pseudo") {
    method_name <- paste0("Augmented Inverse Probability of Treatment",
                          " Weighting: Pseudo-Values")
  } else if (x$method=="tmle") {
    method_name <- "Targeted Maximum Likelihood Estimation"
  } else if (x$method=="tmle_pseudo") {
    method_name <- "Targeted Maximum Likelihood Estimation: Pseudo-Values"
  }

  times_str <- ifelse(is.null(x$call$times), "Event-Specific Times",
                      "User-Supplied Points in Time")

  cat("Confounder Adjusted Cumulative Incidences \n")
  cat("   - Cause of Interest: ", x$call$cause, "\n", sep="")
  cat("   - Method: ", method_name, "\n", sep="")
  cat("   - Times: ", times_str, "\n", sep="")

  if (!is.null(x$boot_data)) {
    cat("   - Bootstrapping: Performed with ", max(x$boot_data$boot),
        " Replications\n", sep="")
  } else {
    cat("   - Bootstrapping: Not Done\n", sep="")
  }

  if (is.null(x$call$conf_int) | !as.logical(as.character(x$call$conf_int))) {
    cat("   - Approximate CI: Not Calculated\n", sep="")
  } else {
    conf_level <- ifelse(is.null(x$call$conf_level), 0.95, x$call$conf_level)
    cat("   - Approximate CI: Calculated with a Confidence level of ",
        conf_level, "\n", sep="")
  }

}
