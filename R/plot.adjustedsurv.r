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

## plot the survival curves
#' @importFrom rlang .data
#' @export
plot.adjustedsurv <- function(x, conf_int=FALSE, max_t=Inf,
                              iso_reg=FALSE, force_bounds=FALSE,
                              use_boot=FALSE, cif=FALSE,
                              color=TRUE, linetype=FALSE, facet=FALSE,
                              line_size=1, line_alpha=1, xlab="Time",
                              ylab="Adjusted Survival Probability",
                              title=NULL, subtitle=NULL, legend.title="Group",
                              legend.position="right",
                              gg_theme=ggplot2::theme_classic(),
                              ylim=NULL, custom_colors=NULL,
                              custom_linetypes=NULL,
                              single_color=NULL, single_linetype=NULL,
                              conf_int_alpha=0.4, steps=TRUE,
                              x_breaks=ggplot2::waiver(), x_n_breaks=NULL,
                              y_breaks=ggplot2::waiver(), y_n_breaks=NULL,
                              additional_layers=list(),
                              median_surv_lines=FALSE, median_surv_size=0.5,
                              median_surv_linetype="dashed",
                              median_surv_color="black", median_surv_alpha=1,
                              median_surv_quantile=0.5,
                              censoring_ind="none",
                              censoring_ind_size=0.5, censoring_ind_alpha=1,
                              censoring_ind_shape=17, censoring_ind_width=NULL,
                              risk_table=FALSE, risk_table_type="n_at_risk",
                              risk_table_stratify=FALSE, risk_table_height=0.25,
                              risk_table_xlab=xlab, risk_table_ylab="default",
                              risk_table_title="default",
                              risk_table_title_size=14,
                              risk_table_title_position="left",
                              risk_table_y_vjust=5, risk_table_theme=gg_theme,
                              risk_table_size=4.2,
                              risk_table_alpha=1, risk_table_color="black",
                              risk_table_family="sans",
                              risk_table_fontface="plain",
                              risk_table_reverse=TRUE,
                              risk_table_stratify_color=TRUE,
                              risk_table_custom_colors=custom_colors,
                              risk_table_use_weights=TRUE,
                              risk_table_digits=1, risk_table_format=TRUE,
                              risk_table_warn=TRUE,
                              risk_table_additional_layers=list(),
                              ...) {
  requireNamespace("ggplot2")

  # get relevant data for the confidence interval
  if (use_boot & is.null(x$boot_adj)) {
    plotdata <- x$adj
  } else if (use_boot) {
    plotdata <- x$boot_adj
  } else {
    plotdata <- x$adj
  }

  # ensure that curves always start at 0
  plotdata <- add_rows_with_zero(plotdata)

  # shortcut to only show curves up to a certain time
  if (is.finite(max_t) && max(plotdata$time) > max_t) {
    max_t_data <- specific_times(plotdata, times=max_t, est="surv",
                                 interpolation=ifelse(steps, "steps", "linear"))
    max_t_data <- max_t_data[!is.na(max_t_data$surv), ]
    plotdata <- plotdata[which(plotdata$time <= max_t), ]
    plotdata <- rbind(plotdata, max_t_data)
  }

  # in some methods estimates can be outside the 0, 1 bounds,
  # if specified set those to 0 or 1 respectively
  if (force_bounds) {
    plotdata <- force_bounds_est(plotdata)
  }

  # apply isotonic regression if specified
  if (iso_reg) {
    plotdata <- iso_reg_est(plotdata)
  }

  # plot CIF instead of survival
  if (cif) {
    plotdata$surv <- 1 - plotdata$surv

    if ("ci_lower" %in% colnames(plotdata)) {
      plotdata$ci_lower <- 1 - plotdata$ci_lower
      plotdata$ci_upper <- 1 - plotdata$ci_upper
    }

    if (ylab=="Adjusted Survival Probability") {
      ylab <- "Adjusted Cumulative Incidence"
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
    line_obj <- ggplot2::geom_step(linewidth=line_size, alpha=line_alpha)
  } else {
    line_obj <- ggplot2::geom_line(linewidth=line_size, alpha=line_alpha)
  }

  # override color using just one color
  if (!is.null(single_color)) {
    if (color) {
      warning("Argument 'color' will be overwritten by argument",
              " 'single_color'.", call.=FALSE)
    }
    line_obj$aes_params$colour <- single_color
  }
  # override color using just one color
  if (!is.null(single_linetype)) {
    if (linetype) {
      warning("Argument 'linetype' will be overwritten by argument",
              " 'single_linetype'.", call.=FALSE)
    }
    line_obj$aes_params$linetype <- single_linetype
  }

  # add steps / lines to plot
  p <- p + line_obj

  # don't use the word "adjusted" with standard Kaplan-Meier
  if (ylab=="Adjusted Survival Probability" & x$method=="km") {
    ylab <- "Survival Probability"
  } else if (ylab=="Adjusted Cumulative Incidence" & x$method=="km") {
    ylab <- "Cumulative Incidence"
  }
  # also warn the user when using steps=FALSE with Kaplan-Meier
  if (x$method=="km" & !steps) {
    warning("Unadjusted Kaplan-Meier estimates should only be drawn as",
            " step functions (steps=TRUE).", call.=FALSE)
  }

  p <- p + gg_theme +
    ggplot2::labs(x=xlab, y=ylab, color=legend.title,
                  linetype=legend.title, fill=legend.title,
                  title=title, subtitle=subtitle) +
    ggplot2::theme(legend.position=legend.position) +
    ggplot2::scale_x_continuous(breaks=x_breaks, n.breaks=x_n_breaks) +
    ggplot2::scale_y_continuous(breaks=y_breaks, n.breaks=y_n_breaks,
                                limits=ylim)

  if (facet) {
    p <- p + ggplot2::facet_wrap(~group)
  }
  if (!is.null(custom_linetypes)) {
    p <- p + ggplot2::scale_linetype_manual(values=custom_linetypes)
  }
  if (!is.null(custom_colors)) {
    p <- p + ggplot2::scale_colour_manual(values=custom_colors)
  }
  if (!is.null(custom_colors)) {
    p <- p + ggplot2::scale_fill_manual(values=custom_colors)
  }

  ## Censoring indicators
  if (censoring_ind!="none") {
    p <- add_censoring_ind(p=p, x=x, plotdata=plotdata, max_t=max_t,
                           steps=steps, color=color, linetype=linetype,
                           ylim=ylim, single_color=single_color,
                           single_linetype=single_linetype,
                           censoring_ind=censoring_ind,
                           censoring_ind_size=censoring_ind_size,
                           censoring_ind_alpha=censoring_ind_alpha,
                           censoring_ind_shape=censoring_ind_shape,
                           censoring_ind_width=censoring_ind_width)
  }

  ## Confidence intervals
  if (conf_int & use_boot & is.null(x$boot_adj)) {
    warning("Cannot use bootstrapped estimates as they were not estimated.",
            " Need bootstrap=TRUE in adjustedsurv() call.", call.=FALSE)
  } else if (conf_int & !use_boot & !"ci_lower" %in% colnames(x$adj)) {
    warning("Cannot draw confidence intervals. Either set 'conf_int=TRUE' in",
            " adjustedsurv() call or use bootstrap estimates.", call.=FALSE)
  } else if (conf_int) {
    p <- add_ci_ribbon(p=p, steps=steps, color=color,
                       single_color=single_color,
                       conf_int_alpha=conf_int_alpha)
  }

  ## Median Survival indicators
  if (median_surv_lines & cif) {
    warning("Cannot draw median survival indicators when using cif=TRUE.",
            call.=FALSE)
  } else if (median_surv_lines) {
    p <- add_median_surv_lines(p=p, x=x, plotdata=plotdata, max_t=max_t,
                               ylim=ylim, steps=steps,
                               median_surv_linetype=median_surv_linetype,
                               median_surv_quantile=median_surv_quantile,
                               median_surv_size=median_surv_size,
                               median_surv_color=median_surv_color,
                               median_surv_alpha=median_surv_alpha)
  }

  # potentially add more stuff to the plot
  if (length(additional_layers) > 0) {
    for (i in seq_len(length(additional_layers))) {
      p <- p + additional_layers[[i]]
    }
  }

  ## adding risk tables
  if (risk_table) {

    check_inputs_risk_table(method=x$method, type=risk_table_type,
                            use_weights=risk_table_use_weights,
                            stratify=risk_table_stratify,
                            warn=risk_table_warn)

    # set correct weights if specified
    if (risk_table_use_weights && !is.null(x$weights) &&
        is.null(x$mids_analyses)) {
      weights <- x$weights
    } else if (risk_table_use_weights && !is.null(x$mids_analyses) &&
               !is.null(x$mids_analyses[[1]]$weights)) {
      weights <- lapply(x$mids_analyses, FUN=function(d){d$weights})
    } else {
      weights <- NULL
    }

    # set correct variable if risk table should be stratified
    if (risk_table_stratify) {
      variable <- x$call$variable
    } else {
      variable <- NULL
    }

    p <- add_risk_table(p_surv=p,
                        data=x$data,
                        event=x$call$event,
                        ev_time=x$call$ev_time,
                        variable=variable,
                        height=risk_table_height,
                        type=risk_table_type,
                        weights=weights,
                        xlab=risk_table_xlab,
                        ylab=risk_table_ylab,
                        title=risk_table_title,
                        title_size=risk_table_title_size,
                        title_position=risk_table_title_position,
                        vjust=risk_table_y_vjust,
                        gg_theme=risk_table_theme,
                        text_size=risk_table_size,
                        text_alpha=risk_table_alpha,
                        text_color=risk_table_color,
                        text_family=risk_table_family,
                        text_fontface=risk_table_fontface,
                        reverse_order=risk_table_reverse,
                        color_groups=risk_table_stratify_color,
                        custom_colors=risk_table_custom_colors,
                        digits=risk_table_digits,
                        text_format=risk_table_format,
                        additional_layers=risk_table_additional_layers)
  }

  return(p)
}

## adding a confidence interval ribbon to the plot
add_ci_ribbon <- function(p, steps, color, single_color, conf_int_alpha) {
  # plot using step-function interpolation
  if (steps) {
    requireNamespace("pammtools")

    ci_map <- ggplot2::aes(ymin=.data$ci_lower,
                           ymax=.data$ci_upper,
                           group=.data$group,
                           fill=.data$group,
                           x=.data$time,
                           y=.data$surv)
    if (!color) {
      ci_map$fill <- NULL
    }

    ribbon <- pammtools::geom_stepribbon(ci_map, alpha=conf_int_alpha,
                                         inherit.aes=FALSE)
    # plot using linear interpolation
  } else {
    ci_map <- ggplot2::aes(ymin=.data$ci_lower,
                           ymax=.data$ci_upper,
                           group=.data$group,
                           fill=.data$group,
                           x=.data$time,
                           y=.data$surv)
    if (!color) {
      ci_map$fill <- NULL
    }

    ribbon <- ggplot2::geom_ribbon(ci_map, alpha=conf_int_alpha,
                                   inherit.aes=FALSE)
  }

  if (!is.null(single_color)) {
    ribbon$aes_params$fill <- single_color
  }
  p <- p + ribbon
  return(p)
}

## adding lines at specific survival time quantile
add_median_surv_lines <- function(p, x, plotdata, max_t, ylim, steps,
                                  median_surv_linetype, median_surv_quantile,
                                  median_surv_size, median_surv_color,
                                  median_surv_alpha) {

  # calculate median survival and add other needed values
  fake_adjsurv <- x
  fake_adjsurv$adj <- plotdata

  median_surv <- adjusted_surv_quantile(fake_adjsurv,
                                        p=median_surv_quantile[[1]],
                                        conf_int=FALSE,
                                        interpolation=ifelse(steps, "steps",
                                                             "linear"))
  median_surv$y <- median_surv_quantile[[1]]
  # set to NA if not in plot
  median_surv$q_surv[median_surv$q_surv > max_t] <- NA

  if (is.null(ylim)) {
    median_surv$yend <- ggplot2::layer_scales(p)$y$range$range[1]
  } else {
    median_surv$yend <- ylim[1]
  }

  # remove if missing
  median_surv <- median_surv[!is.na(median_surv$q_surv), ]

  if (sum(is.na(median_surv$q_surv)) < nrow(median_surv)) {

    median_surv$vert_x <- 0
    median_surv$vert_y <- median_surv_quantile[[1]]
    median_surv$vert_yend <- median_surv_quantile[[1]]

    # draw line on surv_p = 0.5 until it hits the curve
    p <- p + ggplot2::geom_segment(ggplot2::aes(x=.data$vert_x,
                                                xend=.data$q_surv,
                                                y=.data$vert_y,
                                                yend=.data$vert_yend),
                                   inherit.aes=FALSE,
                                   linetype=median_surv_linetype,
                                   linewidth=median_surv_size,
                                   color=median_surv_color,
                                   alpha=median_surv_alpha,
                                   data=median_surv)
    # draw indicator lines from middle to bottom
    p <- p + ggplot2::geom_segment(ggplot2::aes(x=.data$q_surv,
                                                xend=.data$q_surv,
                                                y=median_surv_quantile[[1]],
                                                yend=.data$yend),
                                   inherit.aes=FALSE,
                                   linetype=median_surv_linetype,
                                   linewidth=median_surv_size,
                                   color=median_surv_color,
                                   alpha=median_surv_alpha,
                                   data=median_surv)
  }
  return(p)
}

## create dataset used to plot censoring indicators
get_censoring_ind_data <- function(x, steps, max_t, plotdata) {

  levs <- levels(x$adj$group)

  if (is.null(x$mids_analyses)) {
    data <- x$data
  } else {
    data <- x$data$data
  }

  # keep only relevant data
  data <- data[which(data[, x$call$ev_time] <= max_t), ]

  # create needed data.frame
  cens_dat <- vector(mode="list", length=length(levs))
  for (i in seq_len(length(levs))) {

    # times with censoring
    cens_times <- sort(unique(data[, x$call$ev_time][
      data[, x$call$event]==0 & data[, x$call$variable]==levs[i]]))
    # y axis place to put them
    adjsurv_temp <- plotdata[plotdata$group==levs[i] & !is.na(plotdata$surv), ]

    cens_surv <- read_from_fun(x=cens_times, data=adjsurv_temp,
                               interpolation=ifelse(steps, "steps", "linear"))
    if (length(cens_times)!=0) {
      cens_dat[[i]] <- data.frame(time=cens_times, surv=cens_surv,
                                  group=levs[i])
    }
  }
  cens_dat <- dplyr::bind_rows(cens_dat)
  cens_dat <- cens_dat[!is.na(cens_dat$surv), ]

  return(cens_dat)
}

## add censoring indicators to plot.adjustedsurv()
add_censoring_ind <- function(p, x, plotdata, max_t, steps, color, linetype,
                              single_color, single_linetype, ylim,
                              censoring_ind, censoring_ind_size,
                              censoring_ind_alpha, censoring_ind_shape,
                              censoring_ind_width) {

  if (is.null(censoring_ind_width)) {
    if (is.null(ylim)) {
      ystart <- 1 - ggplot2::layer_scales(p)$y$range$range[1]
    } else {
      ystart <- 1 - ylim[1]
    }
    censoring_ind_width <- ystart * 0.05
  }

  cens_dat <- get_censoring_ind_data(x=x, steps=steps, max_t=max_t,
                                     plotdata=plotdata)

  # either points or lines
  if (censoring_ind=="points") {
    cens_map <- ggplot2::aes(x=.data$time,
                             y=.data$surv,
                             group=.data$group,
                             color=.data$group)
  } else if (censoring_ind=="lines") {
    cens_map <- ggplot2::aes(x=.data$time,
                             y=.data$surv-(censoring_ind_width/2),
                             xend=.data$time,
                             yend=.data$surv+(censoring_ind_width/2),
                             group=.data$group,
                             color=.data$group,
                             linetype=.data$group)
  } else {
    stop("Argument 'censoring_ind' must be either 'none', 'lines' or",
         " 'points'. See documentation.", call.=FALSE)
  }

  if (!color) {
    cens_map$colour <- NULL
  }
  if (!linetype) {
    cens_map$linetype <- NULL
  }

  if (censoring_ind=="points") {
    cens_geom <- ggplot2::geom_point(data=cens_dat, cens_map,
                                     size=censoring_ind_size,
                                     alpha=censoring_ind_alpha,
                                     shape=censoring_ind_shape)
  } else if (censoring_ind=="lines") {
    cens_geom <- ggplot2::geom_segment(data=cens_dat, cens_map,
                                       linewidth=censoring_ind_size,
                                       alpha=censoring_ind_alpha)
  }

  if (!is.null(single_color)) {
    cens_geom$aes_params$colour <- single_color
  }
  if (!is.null(single_linetype)) {
    cens_geom$aes_params$linetype <- single_linetype
  }
  p <- p + cens_geom

  return(p)
}
