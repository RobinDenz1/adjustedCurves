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

## plot the cumulative incidence functions
#' @importFrom rlang .data
#' @export
plot.adjustedcif <- function(x, conf_int=FALSE, max_t=Inf,
                             iso_reg=FALSE, force_bounds=FALSE,
                             use_boot=FALSE, color=TRUE,
                             linetype=FALSE, facet=FALSE,
                             line_size=1, line_alpha=1, xlab="Time",
                             ylab="Adjusted Cumulative Incidence",
                             title=NULL, subtitle=NULL, legend.title="Group",
                             legend.position="right",
                             gg_theme=ggplot2::theme_classic(),
                             ylim=NULL, custom_colors=NULL,
                             custom_linetypes=NULL,
                             single_color=NULL, single_linetype=NULL,
                             conf_int_alpha=0.4, steps=TRUE,
                             censoring_ind="none",
                             censoring_ind_size=0.5, censoring_ind_alpha=1,
                             censoring_ind_shape=17, censoring_ind_width=NULL,
                             ...) {
  requireNamespace("ggplot2")

  # get relevant data for the confidence interval
  if (use_boot & is.null(x$boot_adjcif)) {
    plotdata <- x$adjcif
  } else if (use_boot) {
    plotdata <- x$boot_adjcif
  } else {
    plotdata <- x$adjcif
  }

  # ensure that curves always start at 0
  plotdata <- add_rows_with_zero(plotdata, mode="cif")

  # shortcut to only show curves up to a certain time
  if (is.finite(max_t) && max(plotdata$time) > max_t) {
    max_t_data <- specific_times(plotdata, times=max_t, est="cif",
                                 interpolation=ifelse(steps, "steps", "linear"))
    max_t_data <- max_t_data[!is.na(max_t_data$cif), ]
    plotdata <- plotdata[which(plotdata$time <= max_t), ]
    plotdata <- rbind(plotdata, max_t_data)
  }

  # in some methods estimates can be outside the 0, 1 bounds,
  # if specified set those to 0 or 1 respectively
  if (force_bounds) {
    plotdata <- within(plotdata, {
      cif <- ifelse(cif < 0, 0, cif)
      cif <- ifelse(cif > 1, 1, cif)
    })
  }

  # apply isotonic regression if specified
  if (iso_reg & anyNA(plotdata$cif)) {
    stop("Isotonic Regression cannot be used when there are missing",
         " values in the final CIF estimates.")
  } else if (iso_reg) {
    for (lev in levels(plotdata$group)) {
      temp <- plotdata[plotdata$group==lev, ]

      new <- stats::isoreg(temp$cif)$yf
      plotdata$cif[plotdata$group==lev] <- new

      if (conf_int & "ci_lower" %in% colnames(temp)) {
        diff <- temp$cif - new

        plotdata$ci_lower[plotdata$group==lev] <- temp$ci_lower - diff
        plotdata$ci_upper[plotdata$group==lev] <- temp$ci_upper - diff
      }
    }
  }

  mapping <- ggplot2::aes(x=.data$time, y=.data$cif, color=.data$group,
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

  # don't use the word "adjusted" with standard Aalen-Johansen
  if (ylab=="Adjusted Cumulative Incidence" & x$method=="aalen_johansen") {
    ylab <- "Cumulative Incidence"
  }

  # also warn the user when using steps=FALSE with Aalen-Johansen
  if (x$method=="aalen_johansen" & !steps) {
    warning("Unadjusted Aalen-Johansen estimates should only be drawn as",
            " step functions (steps=TRUE).", call.=FALSE)
  }

  p <- p + gg_theme +
    ggplot2::labs(x=xlab, y=ylab, color=legend.title,
                  linetype=legend.title, fill=legend.title,
                  title=title, subtitle=subtitle) +
    ggplot2::theme(legend.position=legend.position)

  if (facet) {
    p <- p + ggplot2::facet_wrap(~group)
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
  if (censoring_ind!="none") {
    if (is.null(censoring_ind_width)) {
      if (is.null(ylim)) {
        ystart <- 1 - ggplot2::layer_scales(p)$y$range$range[1]
      } else {
        ystart <- 1 - ylim[1]
      }
      censoring_ind_width <- ystart * 0.05
    }

    levs <- levels(plotdata$group)

    # keep only relevant data
    x$data <- x$data[which(x$data[, x$call$ev_time] <= max_t), ]

    # calculate needed data
    cens_dat <- vector(mode="list", length=length(levs))
    for (i in seq_len(length(levs))) {

      x$data <- x$data[which(x$data[, x$call$ev_time] <= max_t), ]
      cens_times <- sort(unique(x$data[, x$call$ev_time][
        x$data[, x$call$event]==0 & x$data[, x$call$variable]==levs[i]]))
      adjcif_temp <- plotdata[plotdata$group==levs[i], ]

      if (steps) {
        read_fun <- read_from_step_function
      } else {
        read_fun <- read_from_linear_function
      }

      cens_cif <- vapply(X=cens_times, FUN=read_fun,
                         FUN.VALUE=numeric(1), data=adjcif_temp,
                         est="cif")
      cens_dat[[i]] <- data.frame(time=cens_times, cif=cens_cif,
                                  group=levs[i])

    }
    cens_dat <- dplyr::bind_rows(cens_dat)
    cens_dat <- cens_dat[!is.na(cens_dat$cif), ]

    # either points or lines
    if (censoring_ind=="points") {
      cens_map <- ggplot2::aes(x=.data$time,
                               y=.data$cif,
                               group=.data$group,
                               color=.data$group)
    } else if (censoring_ind=="lines") {
      cens_map <- ggplot2::aes(x=.data$time,
                               y=.data$cif-(censoring_ind_width/2),
                               xend=.data$time,
                               yend=.data$cif+(censoring_ind_width/2),
                               group=.data$group,
                               color=.data$group,
                               linetype=.data$group)
    } else {
      stop("Argument 'censoring_ind' must be either 'none', 'lines' or",
           " 'points'. See documentation.")
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

  }

  ## Confidence intervals
  if (conf_int & use_boot & is.null(x$boot_adjcif)) {
    warning("Cannot use bootstrapped estimates as they were not estimated.",
            " Need bootstrap=TRUE in adjustedcif() call.")
  } else if (conf_int & !use_boot & !"ci_lower" %in% colnames(plotdata)) {
    warning("Cannot draw confidence intervals. Either set 'conf_int=TRUE' in",
            " adjustedcif() call or use bootstrap estimates.")
  } else if (conf_int) {

    # plot using step-function interpolation
    if (steps) {
      requireNamespace("pammtools")

      ci_map <- ggplot2::aes(ymin=.data$ci_lower,
                             ymax=.data$ci_upper,
                             group=.data$group,
                             fill=.data$group,
                             x=.data$time,
                             y=.data$cif)
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
                             y=.data$cif)
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
  }
  return(p)
}
