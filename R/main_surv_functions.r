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
adjustedsurv <- function(data, variable, ev_time, event, method, conf_int=F,
                         conf_level=0.95, times=NULL, bootstrap=F,
                         n_boot=500, n_cores=1, ...) {

  check_inputs_adjustedsurv(data=data, variable=variable,
                            ev_time=ev_time, event=event, method=method,
                            conf_int=conf_int, conf_level=conf_level,
                            times=times, bootstrap=bootstrap,
                            n_boot=n_boot, ...)

  # define those to remove Notes in devtools::check()
  . <- i <- time <- group <- surv_b <- NULL

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

  levs <- unique(data[,variable])

  # get relevant surv_method function
  surv_fun <- get(paste0("surv_", method))

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
                       "multi_result_class", "adjustedsurv_boot")

      boot_out <- foreach::foreach(i=1:n_boot, .packages=pkgs,
                                   .export=export_objs) %dorng% {

      adjustedsurv_boot(data=data, variable=variable, ev_time=ev_time,
                        event=event, method=method, times_input=times_input,
                        times=times, i=i, surv_fun=surv_fun,
                        levs=levs, ...)
      }
      parallel::stopCluster(cl)

    } else {

      boot_out <- vector(mode="list", length=n_boot)
      for (i in 1:n_boot) {
        boot_out[[i]] <- adjustedsurv_boot(data=data, variable=variable,
                                           ev_time=ev_time, event=event,
                                           method=method, times_input=times_input,
                                           times=times, i=i,
                                           surv_fun=surv_fun, levs=levs, ...)
      }
    }

    # transform into data.frames
    boot_data <- lapply(boot_out, function(x) x$boot_data)
    boot_data_same_t <- lapply(boot_out, function(x) x$boot_data_same_t)

    boot_data <- as.data.frame(dplyr::bind_rows(boot_data))
    boot_data_same_t <- as.data.frame(dplyr::bind_rows(boot_data_same_t))

    # calculate some statistics
    boot_stats <- boot_data_same_t %>%
      dplyr::group_by(., time, group) %>%
      dplyr::summarise(surv=mean(surv_b, na.rm=T),
                       sd=stats::sd(surv_b, na.rm=T),
                       ci_lower=stats::quantile(surv_b,
                                                probs=(1-conf_level)/2,
                                                na.rm=T),
                       ci_upper=stats::quantile(surv_b,
                                                probs=1-((1-conf_level)/2),
                                                na.rm=T),
                       n_boot=sum(!is.na(surv_b)),
                       .groups="drop_last")
  }

  # core of the function
  args <- list(data=data, variable=variable, ev_time=ev_time,
               event=event, conf_int=conf_int, conf_level=conf_level,
               times=times, ...)
  plotdata <- R.utils::doCall(surv_fun, args=args)

  out <- list(adjsurv=plotdata,
              method=method,
              categorical=ifelse(length(unique(data[,variable]))>2, T, F),
              call=match.call())

  if (bootstrap) {
    out$boot_data <- boot_data
    out$boot_data_same_t <- boot_data_same_t
    out$boot_adjsurv <- as.data.frame(boot_stats)
  }

  class(out) <- "adjustedsurv"
  return(out)
}

## perform one bootstrap iteration
adjustedsurv_boot <- function(data, variable, ev_time, event, method,
                              times_input, times, i, surv_fun, levs, ...) {

  indices <- sample(x=rownames(data), size=nrow(data), replace=T)
  boot_samp <- data[indices,]
  # IMPORTANT: keeps SL in tmle methods from failing
  row.names(boot_samp) <- 1:nrow(data)

  # if event specific times are used, use event specific times
  # in bootstrapping as well
  if (is.null(times_input) & method %in% c("km", "iptw_km")) {
    times <- NULL
  } else if (is.null(times_input)) {
    times_boot <- sort(unique(boot_samp[, ev_time][boot_samp[, event]==1]))

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

  if (method %in% c("iptw_km", "iptw_cox", "iptw_pseudo", "aiptw",
                    "aiptw_pseudo")) {
    if (inherits(pass_args$treatment_model, "glm") |
        inherits(pass_args$treatment_model, "multinom")) {
      pass_args$treatment_model <- stats::update(pass_args$treatment_model,
                                                 data=boot_samp, trace=F)
    }
  }

  # call surv_method with correct arguments
  args <- list(data=boot_samp, variable=variable, ev_time=ev_time,
               event=event, conf_int=F, conf_level=0.95, times=times)
  args <- c(args, pass_args)

  adjsurv_boot <- R.utils::doCall(surv_fun, args=args)
  adjsurv_boot$boot <- i

  # read from resulting step function at all t in times
  boot_surv <- vector(mode="list", length=length(levs))
  for (j in 1:length(levs)) {

    if (method %in% c("km", "iptw_km") & is.null(times)) {
      times <- unique(data[,ev_time][data[,variable]==levs[j]])
    }

    surv_boot <- sapply(times, read_from_step_function,
                        step_data=adjsurv_boot[adjsurv_boot$group==levs[j],])

    dat_temp <- data.frame(time=times,
                           surv_b=surv_boot,
                           group=levs[j],
                           boot=i)
    boot_surv[[j]] <- dat_temp
  }
  boot_surv <- as.data.frame(dplyr::bind_rows(boot_surv))

  # output
  result <- multi_result_class()

  result$boot_data <- adjsurv_boot
  result$boot_data_same_t <- boot_surv

  return(result)

}

## plot the survival curves
#' @importFrom rlang .data
#' @export
plot.adjustedsurv <- function(x, draw_ci=F, max_t=Inf,
                              iso_reg=F, force_bounds=F, use_boot=F,
                              color=T, linetype=F, facet=F,
                              line_size=1, xlab="Time",
                              ylab="Adjusted Survival Probability",
                              title=NULL, legend_title="Group",
                              legend_position="right",
                              ylim=NULL, custom_colors=NULL,
                              custom_linetypes=NULL,
                              ci_draw_alpha=0.4, ...) {

  if (!color & !linetype & !facet) {
    stop("Groups must be distinguished with at least one of 'color',",
         "'linetype' or 'facet'. Can't all be FALSE.")
  }

  if (use_boot & is.null(x$boot_adjsurv)) {
    warning("Cannot use bootstrapped estimates as they were not estimated.",
            " Need bootstrap=TRUE in adjustedsurv() call.")
    draw_ci <- F
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
      surv <- ifelse(surv < 0, 0, surv);
      surv <- ifelse(surv > 1, 1, surv);
    })
  }

  # apply isotonic regression if specified
  if (iso_reg) {
    for (lev in levels(plotdata$group)) {
      surv <- plotdata$surv[plotdata$group==lev]
      new <- rev(stats::isoreg(rev(surv))$yf)
      plotdata$surv[plotdata$group==lev] <- new
    }
  }

  mapping <- ggplot2::aes(x=.data$time, y=.data$surv, color=.data$group,
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
    ggplot2::labs(x=xlab, y=ylab, color=legend_title,
                  linetype=legend_title, fill=legend_title) +
    ggplot2::theme(legend.position=legend_position)

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
            " 'adjustedsurv()' call or use bootstrap estimates.")
  }

  if (draw_ci & "ci_lower" %in% colnames(plotdata)) {
    p <- p + pammtools::geom_stepribbon(ggplot2::aes(ymin=.data$ci_lower,
                                            ymax=.data$ci_upper,
                                            fill=.data$group,
                                            x=.data$time,
                                            y=.data$surv),
                                        alpha=ci_draw_alpha, inherit.aes=F)
  }
  return(p)
}
