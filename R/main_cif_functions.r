## Second main function of the package. Is basically a wrapper around
## all cif_method functions, offering additional high level stuff
#' @importFrom dplyr %>%
#' @importFrom doRNG %dorng%
#' @importFrom foreach %dopar%
#' @export
adjustedcif <- function(data, variable, ev_time, event, cause, method,
                        conf_int=F, conf_level=0.95, times=NULL, bootstrap=F,
                        n_boot=500, n_cores=1, ...) {

  check_inputs_adjustedcif(data=data, variable=variable, ev_time=ev_time,
                           event=event, cause=cause, method=method,
                           conf_int=conf_int, conf_level=conf_level,
                           times=times, bootstrap=bootstrap,
                           n_boot=n_boot, ...)

  # define those to remove Notes in devtools::check()
  . <- i <- time <- group <- cif_b <- NULL

  # get event specific times
  times_input <- times
  if (is.null(times)) {
    times <- sort(unique(data[, ev_time][data[, event]>=1]))

    # add zero if not already in there
    if (!0 %in% times) {
      times <- c(0, times)
    }
  }

  levs <- unique(data[,variable])

  # get relevant cif_method function
  cif_fun <- get(paste0("cif_method_", method))

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
                         levs=levs, cause=cause, ...)
      }
      parallel::stopCluster(cl)

    } else {

      boot_out <- vector(mode="list", length=n_boot)
      for (i in 1:n_boot) {
        boot_out[[i]] <- adjustedcif_boot(data=data, variable=variable,
                                          ev_time=ev_time, event=event,
                                          method=method, times_input=times_input,
                                          times=times, i=i, cause=cause,
                                          cif_fun=cif_fun, levs=levs, ...)
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
      dplyr::summarise(cif=mean(cif_b, na.rm=T),
                       sd=stats::sd(cif_b, na.rm=T),
                       ci_lower=stats::quantile(cif_b, probs=1-conf_level,
                                                na.rm=T),
                       ci_upper=stats::quantile(cif_b, probs=conf_level,
                                                na.rm=T),
                       n_boot=sum(!is.na(cif_b)),
                       .groups="drop_last")
  }

  # core of the function
  args <- list(data=data, variable=variable, ev_time=ev_time,
               event=event, conf_int=conf_int, conf_level=conf_level,
               times=times, cause=cause, ...)
  plotdata <- R.utils::doCall(cif_fun, args=args)

  out <- list(adjcif=plotdata,
              method=method,
              categorical=ifelse(length(unique(data[,variable]))>2, T, F),
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
                             times_input, times, i, cif_fun, levs, ...) {

  indices <- sample(x=rownames(data), size=nrow(data), replace=T)
  boot_samp <- data[indices,]

  # if event specific times are used, use event specific times
  # in bootstrapping as well
  if (is.null(times_input)) {
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
                             title=NULL, legend_title="Group",
                             legend_position="right",
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