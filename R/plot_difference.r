
## calculate the curve of the difference
get_diff_curve <- function(x, group_1, group_2, use_boot, conf_level) {

  # silence devtools::check()
  . <- time <- est <- NULL

  # get plotdata out of x
  if (inherits(x, "adjustedsurv")) {
    est_type <- "surv"
    plotdata <- x$adjsurv
  } else if (inherits(x, "adjustedcif")) {
    est_type <- "cif"
    plotdata <- x$adjcif
  }

  # set groups
  if (is.null(group_1) | is.null(group_2)) {
    group_1 <- levels(plotdata$group)[1]
    group_2 <- levels(plotdata$group)[2]
  }

  # keep only data with the groups of interest
  plotdata <- plotdata[plotdata$group==group_1 | plotdata$group==group_2,]
  plotdata$group <- factor(plotdata$group, levels=c(group_1, group_2))

  # observed
  times <- sort(unique(plotdata$time))
  diff_curve <- exact_stepfun_difference(adj=plotdata, times=times,
                                         est=est_type)
  colnames(diff_curve) <- c("time", "est")

  # bootstrapped
  if (use_boot) {

    boot_diff_curves <- vector(mode="list", length=max(x$boot_data$boot))
    for (i in seq_len(max(x$boot_data$boot))) {

      boot_dat <- x$boot_data[x$boot_data$boot==i, ]
      boot_dat$group <- factor(boot_dat$group, levels=c(group_1, group_2))

      boot_diff <- exact_stepfun_difference(adj=boot_dat, times=times,
                                            est=est_type)
      boot_diff$boot <- i
      boot_diff_curves[[i]] <- boot_diff
    }
    boot_diff_curves <- dplyr::bind_rows(boot_diff_curves)
    colnames(boot_diff_curves) <- c("time", "est", "boot")

    # calculate bootstrap confidence intervals
    boot_diff_curve <- boot_diff_curves %>%
      dplyr::group_by(., time) %>%
      dplyr::summarise(se=stats::sd(est, na.rm=TRUE),
                       ci_lower=stats::quantile(est,
                                                probs=(1-conf_level)/2,
                                                na.rm=TRUE),
                       ci_upper=stats::quantile(est,
                                                probs=1-((1-conf_level)/2),
                                                na.rm=TRUE),
                       n_boot=sum(!is.na(est)),
                       .groups="drop_last")
    boot_diff_curve$est <- diff_curve$est
    diff_curve <- boot_diff_curve
  }

  output <- list(diff_curve=diff_curve,
                 group_1=group_1,
                 group_2=group_2,
                 est_type=est_type)

  return(output)
}

## plot the difference between two adjusted survival curves
#' @importFrom rlang .data
#' @export
plot_difference <- function(x, group_1=NULL, group_2=NULL, conf_int=FALSE,
                            conf_level=0.95, type="steps", max_t=Inf,
                            size=0.7, color="black", linetype="solid", alpha=1,
                            conf_int_alpha=0.4, points_ci_size=NULL,
                            points_ci_width=NULL, xlab="Time", ylab=NULL,
                            title=NULL, subtitle=NULL,
                            gg_theme=ggplot2::theme_classic(),
                            line_at_0=TRUE, line_at_0_size=0.7,
                            line_at_0_color="grey", line_at_0_linetype="dashed",
                            line_at_0_alpha=1,
                            loess_smoother=FALSE, loess_span=0.75,
                            loess_color=color, loess_size=size,
                            loess_linetype="dashed", loess_alpha=alpha,
                            ...) {
  requireNamespace("ggplot2")

  check_inputs_plot_difference(x=x, group_1=group_1, group_2=group_2,
                               conf_int=conf_int, type=type, max_t=max_t)

  # get relevant data
  diff_obj <- get_diff_curve(x=x, group_1=group_1, group_2=group_2,
                             use_boot=conf_int, conf_level=conf_level)
  plotdata <- diff_obj$diff_curve
  plotdata <- plotdata[which(plotdata$time <= max_t), ]

  # initialize plot
  p <- ggplot2::ggplot(plotdata, ggplot2::aes(x=.data$time, y=.data$est))

  # add line at 0 if specified
  if (line_at_0) {
    p <- p + ggplot2::geom_hline(yintercept=0, size=line_at_0_size,
                                 color=line_at_0_color,
                                 linetype=line_at_0_linetype,
                                 alpha=line_at_0_alpha)
  }

  # plot difference as step function
  if (type=="steps") {
    if (conf_int) {
      requireNamespace("pammtools")

      p <- p + pammtools::geom_stepribbon(ggplot2::aes(ymin=.data$ci_lower,
                                                       ymax=.data$ci_upper,
                                                       x=.data$time,
                                                       y=.data$est),
                                          alpha=conf_int_alpha,
                                          fill=color,
                                          inherit.aes=FALSE)
    }
    p <- p + ggplot2::geom_step(size=size, color=color, linetype=linetype,
                                alpha=alpha)
  # plot difference using linear interpolation
  } else if (type=="lines") {
    p <- p + ggplot2::geom_line(size=size, color=color, linetype=linetype,
                                alpha=alpha)
    if (conf_int) {
      p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin=.data$ci_lower,
                                                 ymax=.data$ci_upper,
                                                 x=.data$time,
                                                 y=.data$est),
                                    alpha=conf_int_alpha,
                                    fill=color,
                                    inherit.aes=FALSE)
    }
  } else if (type=="points") {
    if (conf_int) {

      if (is.null(points_ci_size)) {
        points_ci_size <- max(plotdata$time / 100)
      }

      if (is.null(points_ci_width)) {
        points_ci_width <- max(plotdata$time / 100)
      }

      p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin=.data$ci_lower,
                                                   ymax=.data$ci_upper),
                                      size=points_ci_size,
                                      width=points_ci_width,
                                      color=color)
    }
    p <- p + ggplot2::geom_point(size=size, color=color, alpha=alpha)
  }

  # add loess smoother line
  if (loess_smoother) {
    p <- p + ggplot2::geom_line(stat="smooth", method="loess",
                                formula=y ~ x, se=FALSE,
                                color=loess_color, size=loess_size,
                                alpha=loess_alpha, linetype=loess_linetype,
                                span=loess_span)
  }

  # generate default label
  if (is.null(group_1) | is.null(group_2)) {
    group_1 <- diff_obj$group_1
    group_2 <- diff_obj$group_2
  }

  if (is.null(ylab) & diff_obj$est_type=="surv") {
    ylab <- bquote(hat(S)[.(group_1)](t) - hat(S)[.(group_2)](t))
  } else if (is.null(ylab) & diff_obj$est_type=="cif") {
    ylab <- bquote(hat(F)[.(group_1)](t) - hat(F)[.(group_2)](t))
  }

  p <- p + gg_theme
  p <- p + ggplot2::labs(x=xlab, y=ylab, title=title, subtitle=subtitle)

  return(p)
}
