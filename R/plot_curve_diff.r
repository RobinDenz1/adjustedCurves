
## plot the difference between two adjusted survival curves
#' @importFrom rlang .data
#' @export
plot_curve_diff <- function(x, group_1=NULL, group_2=NULL, conf_int=FALSE,
                            conf_level=0.95, type="steps", times=NULL,
                            max_t=Inf, use_boot=FALSE, size=0.7, color="black",
                            linetype="solid", alpha=1,
                            conf_int_alpha=0.4, points_ci_size=NULL,
                            points_ci_width=NULL, xlab="Time", ylab=NULL,
                            title=NULL, subtitle=NULL,
                            gg_theme=ggplot2::theme_classic(),
                            line_at_ref=TRUE, line_at_ref_size=0.7,
                            line_at_ref_color="grey",
                            line_at_ref_linetype="dashed",
                            line_at_ref_alpha=1,
                            loess_smoother=FALSE, loess_span=0.75,
                            loess_color=color, loess_size=size,
                            loess_linetype="dashed", loess_alpha=alpha,
                            test=NULL, integral_from=0, integral_to=NULL,
                            p_value=FALSE, integral=FALSE,
                            interval=FALSE, text_pos_x="left",
                            text_pos_y="bottom", text_size=3.5,
                            text_family="serif", text_fontface="italic",
                            text_color="black", text_alpha=1,
                            text_digits=3, text_format_p=TRUE,
                            fill_area=FALSE, area_color="blue", area_alpha=0.4,
                            fill_only_interval=TRUE,
                            ...) {
  requireNamespace("ggplot2")

  check_inputs_plot_difference(x=x, group_1=group_1, group_2=group_2,
                               conf_int=conf_int, type=type, max_t=max_t,
                               test=test, integral_from=integral_from,
                               integral_to=integral_to, p_value=p_value,
                               integral=integral, use_boot=use_boot)

  # object specific stuff
  adj_data <- x$adj
  if (inherits(x, "adjustedsurv")) {
    mode <- "surv"
  } else {
    mode <- "cif"
  }

  # what kind of interpolation to use
  if (type=="lines") {
    interpolation <- "linear"
  } else {
    interpolation <- "steps"
  }

  # get relevant data
  plotdata <- adjusted_curve_diff(adj=x, group_1=group_1, group_2=group_2,
                                  conf_int=conf_int, conf_level=conf_level,
                                  interpolation=interpolation, times=times,
                                  use_boot=use_boot)
  plotdata <- plotdata[which(!is.na(plotdata$diff)), ]
  plotdata <- plotdata[which(plotdata$time <= max_t), ]

  # initialize plot
  p <- ggplot2::ggplot(plotdata, ggplot2::aes(x=.data$time, y=.data$diff))

  # add line at 0 if specified
  if (line_at_ref) {
    p <- p + ggplot2::geom_hline(yintercept=0,
                                 linewidth=line_at_ref_size,
                                 color=line_at_ref_color,
                                 linetype=line_at_ref_linetype,
                                 alpha=line_at_ref_alpha)
  }

  # plot difference as step function
  if (type=="steps") {
    if (conf_int) {
      requireNamespace("pammtools")

      p <- p + pammtools::geom_stepribbon(ggplot2::aes(ymin=.data$ci_lower,
                                                       ymax=.data$ci_upper,
                                                       x=.data$time,
                                                       y=.data$diff),
                                          alpha=conf_int_alpha,
                                          fill=color,
                                          inherit.aes=FALSE)
    }
    p <- p + ggplot2::geom_step(linewidth=size, color=color, linetype=linetype,
                                alpha=alpha)
  # plot difference using linear interpolation
  } else if (type=="lines") {
    p <- p + ggplot2::geom_line(linewidth=size, color=color, linetype=linetype,
                                alpha=alpha)
    if (conf_int) {
      p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin=.data$ci_lower,
                                                 ymax=.data$ci_upper,
                                                 x=.data$time,
                                                 y=.data$diff),
                                    alpha=conf_int_alpha,
                                    fill=color,
                                    inherit.aes=FALSE)
    }
  # plot difference using points and maybe errorbars
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
                                      linewidth=points_ci_size,
                                      width=points_ci_width,
                                      color=color)
    }
    p <- p + ggplot2::geom_point(size=size, color=color, alpha=alpha)
  }

  # add loess smoother line
  if (loess_smoother) {
    p <- p + ggplot2::geom_line(stat="smooth", method="loess",
                                formula=y ~ x, se=FALSE,
                                color=loess_color, linewidth=loess_size,
                                alpha=loess_alpha, linetype=loess_linetype,
                                span=loess_span)
  }

  # generate default label
  if (is.null(group_1) | is.null(group_2)) {
    group_1 <- levels(adj_data$group)[1]
    group_2 <- levels(adj_data$group)[2]
  }

  if (is.null(ylab) & mode=="surv") {
    ylab <- bquote(hat(S)[.(group_1)](t) - hat(S)[.(group_2)](t))
  } else if (is.null(ylab) & mode=="cif") {
    ylab <- bquote(hat(F)[.(group_1)](t) - hat(F)[.(group_2)](t))
  }

  p <- p + gg_theme
  p <- p + ggplot2::labs(x=xlab, y=ylab, title=title, subtitle=subtitle)

  # perform test here if "test" is NULL and p-val is wanted
  if (p_value & is.null(test)) {
    test <- adjusted_curve_test(adj=x, to=integral_to, from=integral_from,
                                conf_level=conf_level,
                                interpolation=interpolation,
                                group_1=group_1, group_2=group_2)
    p_val <- test$p_value
  # if only the integral should be printed, calculate this only
  } else if (integral & is.null(test)) {
    area <- exact_integral(data=plotdata, from=integral_from,
                           to=integral_to, interpolation=interpolation,
                           est="diff")
  }

  # add p-value and other text to the plot, based on a test performed by
  # the adjusted_curve_test function
  if (p_value | integral) {
    requireNamespace("ggpp")

    if (!is.null(test)) {
      p_val <- test$p_value
      area <- test$observed_diff_integral
      to <- test$call$to

      if (!is.numeric(test$call$from)) {
        from <- 0
      } else {
        from <- test$call$from
      }
    } else {
      to <- integral_to
      from <- integral_from
    }

    # create label
    lab <- ""
    if (p_value & text_format_p) {
      p_form <- format.pval(p_val, digits=text_digits, eps=0.01)
      if (startsWith(p_form, "<")) {
        lab <- paste0(lab, "p ", p_form)
      } else {
        lab <- paste0(lab, "p = ", p_form)
      }
    } else if (p_value) {
      lab <- paste0(lab, "p = ", round(p_val, text_digits))
    }
    if (integral) {
      lab <- paste0(lab, "\nArea = ", round(area, text_digits))
    }
    if (interval) {
      lab <- paste0(lab, "\nInterval: [", from, ", ", to, "]")
    }

    # put together
    p_dat <- data.frame(x=text_pos_x,
                        y=text_pos_y,
                        label=lab)
    p <- p + ggpp::geom_text_npc(data=p_dat,
                                 ggplot2::aes(npcx=.data$x, npcy=.data$y,
                                              label=.data$label),
                                 size=text_size,
                                 family=text_family,
                                 fontface=text_fontface,
                                 color=text_color,
                                 alpha=text_alpha)
  }

  # restrict area to interval used for calculations
  if (fill_area & fill_only_interval) {
    restricted_times <- sort(unique(plotdata$time))
    restricted_times <- restricted_times[restricted_times <= to &
                                         restricted_times >= from]

    restricted_est <- read_from_fun(x=restricted_times, est="diff",
                                    interpolation=interpolation, data=plotdata)
    plotdata_r <- data.frame(time=restricted_times,
                             diff=restricted_est)

  } else {
    plotdata_r <- plotdata
  }

  # add color to non-zero area if specified
  if (type=="lines" & fill_area) {
    p <- p + ggplot2::geom_area(data=plotdata_r, fill=area_color,
                                alpha=area_alpha)
  } else if (type=="steps" & fill_area) {
    p <- p + pammtools::geom_stepribbon(ggplot2::aes(ymin=0,
                                                     ymax=.data$diff,
                                                     x=.data$time,
                                                     y=.data$diff),
                                        fill=area_color,
                                        alpha=area_alpha,
                                        data=plotdata_r)
  } else if ((type=="none" | type=="points") & fill_area) {
    warning("'fill_area' can only be used with type='lines' and",
            " type='steps'.")
  }

  return(p)
}

## plot the ratio of two adjusted survival curves
#' @importFrom rlang .data
#' @export
plot_curve_ratio <- function(x, group_1=NULL, group_2=NULL, conf_int=FALSE,
                             conf_level=0.95, type="steps", times=NULL,
                             max_t=Inf, use_boot=FALSE, size=0.7, color="black",
                             linetype="solid", alpha=1,
                             conf_int_alpha=0.4, xlab="Time", ylab=NULL,
                             title=NULL, subtitle=NULL,
                             gg_theme=ggplot2::theme_classic(),
                             line_at_ref=TRUE, line_at_ref_size=0.7,
                             line_at_ref_color="grey",
                             line_at_ref_linetype="dashed",
                             line_at_ref_alpha=1, ...) {
  requireNamespace("ggplot2")

  check_inputs_plot_difference(x=x, group_1=group_1, group_2=group_2,
                               conf_int=conf_int, type=type, max_t=max_t,
                               test=NULL, integral_from=0,
                               integral_to=NULL, p_value=FALSE,
                               integral=FALSE, use_boot=use_boot)

  # object specific stuff
  adj_data <- x$adj
  if (inherits(x, "adjustedsurv")) {
    mode <- "surv"
  } else {
    mode <- "cif"
  }

  # what kind of interpolation to use
  if (type=="lines") {
    interpolation <- "linear"
  } else {
    interpolation <- "steps"
  }

  # get relevant data
  plotdata <- adjusted_curve_ratio(adj=x, group_1=group_1, group_2=group_2,
                                   conf_int=conf_int, conf_level=conf_level,
                                   interpolation=interpolation, times=times,
                                   use_boot=use_boot)
  plotdata <- plotdata[which(!is.na(plotdata$ratio)), ]
  plotdata <- plotdata[which(plotdata$time <= max_t), ]

  # initialize plot
  p <- ggplot2::ggplot(plotdata, ggplot2::aes(x=.data$time, y=.data$ratio))

  # add line at 0 if specified
  if (line_at_ref) {
    p <- p + ggplot2::geom_hline(yintercept=1,
                                 linewidth=line_at_ref_size,
                                 color=line_at_ref_color,
                                 linetype=line_at_ref_linetype,
                                 alpha=line_at_ref_alpha)
  }

  # plot difference as step function
  if (type=="steps") {
    if (conf_int) {
      requireNamespace("pammtools")

      p <- p + pammtools::geom_stepribbon(ggplot2::aes(ymin=.data$ci_lower,
                                                       ymax=.data$ci_upper,
                                                       x=.data$time,
                                                       y=.data$ratio),
                                          alpha=conf_int_alpha,
                                          fill=color,
                                          inherit.aes=FALSE)
    }
    p <- p + ggplot2::geom_step(linewidth=size, color=color, linetype=linetype,
                                alpha=alpha)
    # plot difference using linear interpolation
  } else if (type=="lines") {
    p <- p + ggplot2::geom_line(linewidth=size, color=color, linetype=linetype,
                                alpha=alpha)
    if (conf_int) {
      p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin=.data$ci_lower,
                                                 ymax=.data$ci_upper,
                                                 x=.data$time,
                                                 y=.data$ratio),
                                    alpha=conf_int_alpha,
                                    fill=color,
                                    inherit.aes=FALSE)
    }
  }

  # generate default label
  if (is.null(group_1) | is.null(group_2)) {
    group_1 <- levels(adj_data$group)[1]
    group_2 <- levels(adj_data$group)[2]
  }

  if (is.null(ylab) & mode=="surv") {
    ylab <- bquote(hat(S)[.(group_1)](t) / hat(S)[.(group_2)](t))
  } else if (is.null(ylab) & mode=="cif") {
    ylab <- bquote(hat(F)[.(group_1)](t) / hat(F)[.(group_2)](t))
  }

  p <- p + gg_theme
  p <- p + ggplot2::labs(x=xlab, y=ylab, title=title, subtitle=subtitle)

  return(p)
}
