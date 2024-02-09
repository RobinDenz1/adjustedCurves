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

## plot adjusted restricted mean survival time curve
#' @export
plot_rmst_curve <- function(adjsurv, times=NULL, conf_int=FALSE,
                            conf_level=0.95, interpolation="steps", max_t=Inf,
                            color=TRUE, linetype=FALSE, facet=FALSE,
                            size=1, alpha=1, xlab="Time", ylab="RMST",
                            title=NULL, subtitle=NULL, legend.title="Group",
                            legend.position="right",
                            gg_theme=ggplot2::theme_classic(),
                            custom_colors=NULL, custom_linetypes=NULL,
                            conf_int_alpha=0.4, ...) {

  if (!inherits(adjsurv, "adjustedsurv")) {
    stop("'adjsurv' must be an adjustedsurv object created using the",
         " adjustedsurv function.")
  }

  plot_auc_curve(estimate="rmst", adj=adjsurv, times=times, conf_int=conf_int,
                 interpolation=interpolation, max_t=max_t, color=color,
                 linetype=linetype, facet=facet, size=size, alpha=alpha,
                 xlab=xlab, ylab=ylab, title=title,
                 subtitle=subtitle, legend.title=legend.title,
                 legend.position=legend.position, gg_theme=gg_theme,
                 custom_colors=custom_colors, custom_linetypes=custom_linetypes,
                 conf_int_alpha=conf_int_alpha, conf_level=conf_level, ...)
}

## plot adjusted restricted mean time lost curve
#' @export
plot_rmtl_curve <- function(adj, times=NULL, conf_int=FALSE,
                            conf_level=0.95, interpolation="steps", max_t=Inf,
                            color=TRUE, linetype=FALSE, facet=FALSE,
                            size=1, alpha=1, xlab="Time", ylab="RMTL",
                            title=NULL, subtitle=NULL, legend.title="Group",
                            legend.position="right",
                            gg_theme=ggplot2::theme_classic(),
                            custom_colors=NULL, custom_linetypes=NULL,
                            conf_int_alpha=0.4, ...) {

  if (!inherits(adj, c("adjustedsurv", "adjustedcif"))) {
    stop("'adj' must be either an adjustedsurv object created using the",
         " adjustedsurv function or an adjustedcif object created using",
         " the adjustedcif function.")
  }

  plot_auc_curve(estimate="rmtl", adj=adj, times=times, conf_int=conf_int,
                 interpolation=interpolation, max_t=max_t, color=color,
                 linetype=linetype, facet=facet, size=size, alpha=alpha,
                 xlab=xlab, ylab=ylab, title=title,
                 subtitle=subtitle, legend.title=legend.title,
                 legend.position=legend.position, gg_theme=gg_theme,
                 custom_colors=custom_colors, custom_linetypes=custom_linetypes,
                 conf_int_alpha=conf_int_alpha, conf_level=conf_level, ...)
}
