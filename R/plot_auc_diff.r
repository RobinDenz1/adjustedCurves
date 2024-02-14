

plot_auc_diff <- function(adj, estimate, times=NULL, conf_int=FALSE,
                          conf_level=0.95, interpolation="steps", max_t=Inf,
                          type, group_1=NULL, group_2=NULL,
                          color="black", linetype="solid",
                          size=1, alpha=1, xlab="Time", ylab,
                          title=NULL, subtitle=NULL,
                          gg_theme=ggplot2::theme_classic(),
                          conf_int_alpha=0.4,
                          line_at_ref=TRUE, line_at_ref_size=0.7,
                          line_at_ref_color="grey",
                          line_at_ref_linetype="dashed",
                          line_at_ref_alpha=1, ...) {
  # default times
  if (is.null(times)) {
    if (inherits(adj, "adjustedcif")) {
      times <- sort(unique(adj$adjcif$time))
    } else if (inherits(adj, "adjustedsurv")) {
      times <- sort(unique(adj$adjsurv$time))
    }
    times <- times[times > 0]
  }
  if (is.finite(max_t)) {
    times <- times[times < max_t]
    times <- c(times, max_t)
  }

  if (type=="diff") {
    difference <- TRUE
    ratio <- FALSE
    ref_value <- 0
  } else if (type=="ratio") {
    difference <- FALSE
    ratio <- TRUE
    ref_value <- 1
  }

  # calculate rmst curve
  if (estimate=="rmst") {
    plotdata <- adjusted_rmst(adjsurv=adj, from=0, to=times, conf_int=conf_int,
                              conf_level=conf_level, group_1=group_1,
                              group_2=group_2, interpolation=interpolation,
                              ratio=ratio, difference=difference)
  } else if (estimate=="rmtl") {
    plotdata <- adjusted_rmtl(adj=adj, from=0, to=times, conf_int=conf_int,
                              conf_level=conf_level, group_1=group_1,
                              group_2=group_2, interpolation=interpolation,
                              ratio=ratio, difference=difference)
  }
  plotdata$se <- NULL

  # get one consistent name
  if (conf_int) {
    colnames(plotdata) <- c("to", "est", "ci_lower", "ci_upper", "p_value")
  } else {
    colnames(plotdata) <- c("to", "est")
  }

  # remove NAs
  plotdata <- plotdata[!is.na(plotdata$est), ]

  # start plotting
  p <- ggplot2::ggplot(plotdata, ggplot2::aes(x=.data$to, y=.data$est))

  # add line at 0 if specified
  if (line_at_ref) {
    p <- p + ggplot2::geom_hline(yintercept=ref_value,
                                 linewidth=line_at_ref_size,
                                 color=line_at_ref_color,
                                 linetype=line_at_ref_linetype,
                                 alpha=line_at_ref_alpha)
  }

  p <- p +
    ggplot2::geom_line(linewidth=size, alpha=alpha) +
    ggplot2::labs(x=xlab, y=ylab, title=title, subtitle=subtitle) +
    gg_theme

  if (conf_int) {
    ci_map <- ggplot2::aes(ymin=.data$ci_lower,
                           ymax=.data$ci_upper,
                           x=.data$to,
                           y=.data$est)
    p <- p + ggplot2::geom_ribbon(ci_map,
                                  alpha=conf_int_alpha,
                                  inherit.aes=FALSE)
  }
  return(p)
}

