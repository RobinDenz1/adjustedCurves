
## generalized plot function for both rmst and rmtl, called internally
plot_auc_curve <- function(adj, estimate, times=NULL, conf_int=FALSE,
                           conf_level=0.95, interpolation="steps", max_t=Inf,
                           color=TRUE, linetype=FALSE, facet=FALSE,
                           size=1, alpha=1, xlab="Time", ylab,
                           title=NULL, subtitle=NULL, legend.title="Group",
                           legend.position="right",
                           gg_theme=ggplot2::theme_classic(),
                           custom_colors=NULL, custom_linetypes=NULL,
                           conf_int_alpha=0.4, ...) {

  check_inputs_auc_curve(times=times, max_t=max_t, color=color,
                         linetype=linetype, facet=facet)

  # default times
  if (is.null(times)) {
    times <- sort(unique(adj$adj$time))
    times <- times[times > 0]
  }
  if (is.finite(max_t)) {
    times <- times[times < max_t]
    times <- c(times, max_t)
  }

  # calculate rmst curve
  if (estimate=="rmst") {
    plotdata <- adjusted_rmst(adjsurv=adj, from=0, to=times, conf_int=conf_int,
                              conf_level=conf_level,
                              interpolation=interpolation)
  } else if (estimate=="rmtl") {
    plotdata <- adjusted_rmtl(adj=adj, from=0, to=times, conf_int=conf_int,
                              conf_level=conf_level,
                              interpolation=interpolation)
  }
  plotdata <- dplyr::bind_rows(plotdata)
  plotdata$se <- NULL

  # get one consistent name
  if (conf_int) {
    colnames(plotdata) <- c("to", "group", "auc", "ci_lower", "ci_upper",
                            "n_boot")
  } else {
    colnames(plotdata) <- c("to", "group", "auc")
  }

  # remove NAs
  plotdata <- plotdata[!is.na(plotdata$auc), ]

  # start plotting
  mapping <- ggplot2::aes(x=.data$to, y=.data$auc, color=.data$group,
                          linetype=.data$group, group=.data$group)
  if (!linetype) {
    mapping$linetype <- NULL
  }
  if (!color) {
    mapping$colour <- NULL
  }

  p <- ggplot2::ggplot(plotdata, mapping) +
    ggplot2::geom_line(linewidth=size, alpha=alpha) +
    ggplot2::labs(x=xlab, y=ylab, title=title, subtitle=subtitle,
                  color=legend.title, linetype=legend.title,
                  fill=legend.title) +
    gg_theme +
    ggplot2::theme(legend.position=legend.position)

  if (facet) {
    p <- p + ggplot2::facet_wrap(~group)
  }
  if (!is.null(custom_colors)) {
    p <- p + ggplot2::scale_colour_manual(values=custom_colors)
    p <- p + ggplot2::scale_fill_manual(values=custom_colors)
  }
  if (!is.null(custom_linetypes)) {
    p <- p + ggplot2::scale_linetype_manual(values=custom_linetypes)
  }
  if (conf_int) {
    ci_map <- ggplot2::aes(ymin=.data$ci_lower,
                           ymax=.data$ci_upper,
                           group=.data$group,
                           fill=.data$group,
                           x=.data$to,
                           y=.data$auc)
    if (!color) {
      ci_map$fill <- NULL
    }
    p <- p + ggplot2::geom_ribbon(ci_map, alpha=conf_int_alpha,
                                  inherit.aes=FALSE)
  }
  return(p)
}

## internal plot function called when the difference of ratio
## of two RMST / RMTL curves should be plotted instead
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
    times <- sort(unique(adj$adj$time))
    times <- times[times > 0]
  }
  if (is.finite(max_t)) {
    times <- times[times < max_t]
    times <- c(times, max_t)
  }

  if (type=="diff") {
    ref_value <- 0
  } else if (type=="ratio") {
    ref_value <- 1
  }

  # calculate rmst curve
  if (estimate=="rmst") {
    plotdata <- adjusted_rmst(adjsurv=adj, from=0, to=times, conf_int=conf_int,
                              conf_level=conf_level, group_1=group_1,
                              group_2=group_2, interpolation=interpolation,
                              contrast=type)
  } else if (estimate=="rmtl") {
    plotdata <- adjusted_rmtl(adj=adj, from=0, to=times, conf_int=conf_int,
                              conf_level=conf_level, group_1=group_1,
                              group_2=group_2, interpolation=interpolation,
                              contrast=type)
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
    ggplot2::geom_line(linewidth=size, alpha=alpha, color=color,
                       linetype=linetype) +
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

## plot adjusted restricted mean survival time curve
#' @export
plot_rmst_curve <- function(adjsurv, times=NULL, conf_int=FALSE,
                            conf_level=0.95, interpolation="steps",
                            contrast="none", group_1=NULL, group_2=NULL,
                            max_t=Inf, color=TRUE, linetype=FALSE, facet=FALSE,
                            size=1, alpha=1, xlab="Time", ylab="RMST",
                            title=NULL, subtitle=NULL, legend.title="Group",
                            legend.position="right",
                            gg_theme=ggplot2::theme_classic(),
                            custom_colors=NULL, custom_linetypes=NULL,
                            conf_int_alpha=0.4,
                            line_at_ref=TRUE, line_at_ref_size=0.7,
                            line_at_ref_color="grey",
                            line_at_ref_linetype="dashed",
                            line_at_ref_alpha=1, ...) {

  if (!inherits(adjsurv, "adjustedsurv")) {
    stop("'adjsurv' must be an adjustedsurv object created using the",
         " adjustedsurv function.")
  }

  if (!contrast %in% c("diff", "ratio")) {
    p <- plot_auc_curve(estimate="rmst", adj=adjsurv, times=times,
                        conf_int=conf_int, interpolation=interpolation,
                        max_t=max_t, color=color,
                        linetype=linetype, facet=facet, size=size,
                        alpha=alpha, xlab=xlab, ylab=ylab, title=title,
                        subtitle=subtitle, legend.title=legend.title,
                        legend.position=legend.position, gg_theme=gg_theme,
                        custom_colors=custom_colors,
                        custom_linetypes=custom_linetypes,
                        conf_int_alpha=conf_int_alpha, conf_level=conf_level,
                        ...)
  } else {

    if (!is.character(color)) {
      color <- "black"
    }
    if (!is.character(linetype)) {
      linetype <- "solid"
    }

    check_inputs_auc_diff(times=times, max_t=max_t, color=color,
                          contrast=contrast, linetype=linetype,
                          line_at_ref=line_at_ref)

    p <- plot_auc_diff(adj=adjsurv, estimate="rmst", times=times,
                       conf_int=conf_int, conf_level=conf_level,
                       interpolation=interpolation, max_t=max_t,
                       type=contrast, group_1=group_1, group_2=group_2,
                       color=color, linetype=linetype,
                       size=size, alpha=alpha, xlab=xlab, ylab=ylab,
                       title=title, subtitle=subtitle,
                       gg_theme=gg_theme, conf_int_alpha=conf_int_alpha,
                       line_at_ref=line_at_ref,
                       line_at_ref_size=line_at_ref_size,
                       line_at_ref_color=line_at_ref_color,
                       line_at_ref_linetype=line_at_ref_linetype,
                       line_at_ref_alpha=line_at_ref_alpha, ...)
  }
  return(p)
}

## plot adjusted restricted mean time lost curve
#' @export
plot_rmtl_curve <- function(adj, times=NULL, conf_int=FALSE,
                            conf_level=0.95, interpolation="steps",
                            contrast="none", group_1=NULL, group_2=NULL,
                            max_t=Inf, color=TRUE, linetype=FALSE, facet=FALSE,
                            size=1, alpha=1, xlab="Time", ylab="RMTL",
                            title=NULL, subtitle=NULL, legend.title="Group",
                            legend.position="right",
                            gg_theme=ggplot2::theme_classic(),
                            custom_colors=NULL, custom_linetypes=NULL,
                            conf_int_alpha=0.4,
                            line_at_ref=TRUE, line_at_ref_size=0.7,
                            line_at_ref_color="grey",
                            line_at_ref_linetype="dashed",
                            line_at_ref_alpha=1, ...) {

  if (!inherits(adj, c("adjustedsurv", "adjustedcif"))) {
    stop("'adj' must be either an adjustedsurv object created using the",
         " adjustedsurv function or an adjustedcif object created using",
         " the adjustedcif function.")
  }

  if (!contrast %in% c("diff", "ratio")) {
    p <- plot_auc_curve(estimate="rmtl", adj=adj, times=times,
                        conf_int=conf_int, interpolation=interpolation,
                        max_t=max_t, color=color, linetype=linetype,
                        facet=facet, size=size, alpha=alpha,
                        xlab=xlab, ylab=ylab, title=title,
                        subtitle=subtitle, legend.title=legend.title,
                        legend.position=legend.position, gg_theme=gg_theme,
                        custom_colors=custom_colors,
                        custom_linetypes=custom_linetypes,
                        conf_int_alpha=conf_int_alpha,
                        conf_level=conf_level, ...)
  } else {

    if (!is.character(color)) {
      color <- "black"
    }
    if (!is.character(linetype)) {
      linetype <- "solid"
    }

    check_inputs_auc_diff(times=times, max_t=max_t, color=color,
                          contrast=contrast, linetype=linetype,
                          line_at_ref=line_at_ref)

    p <- plot_auc_diff(adj=adj, estimate="rmtl", times=times,
                       conf_int=conf_int, conf_level=conf_level,
                       interpolation=interpolation, max_t=max_t,
                       type=contrast, group_1=group_1, group_2=group_2,
                       color=color, linetype=linetype,
                       size=size, alpha=alpha, xlab=xlab, ylab=ylab,
                       title=title, subtitle=subtitle,
                       gg_theme=gg_theme, conf_int_alpha=conf_int_alpha,
                       line_at_ref=line_at_ref,
                       line_at_ref_size=line_at_ref_size,
                       line_at_ref_color=line_at_ref_color,
                       line_at_ref_linetype=line_at_ref_linetype,
                       line_at_ref_alpha=line_at_ref_alpha, ...)
  }
  return(p)
}
