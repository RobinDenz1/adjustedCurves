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

## calculate the curve of the difference including
## confidence intervals and p-values
#' @export
adjusted_curve_diff <- function(adj, group_1=NULL, group_2=NULL, times=NULL,
                                conf_int=FALSE, conf_level=0.95,
                                use_boot=FALSE, interpolation="steps") {

  check_inputs_adj_diff(adj=adj, group_1=group_1, group_2=group_2,
                        conf_int=conf_int, use_boot=use_boot)

  # silence devtools::check()
  . <- time <- est <- NULL

  # get plotdata out of adj
  if (inherits(adj, "adjustedsurv")) {
    est_type <- "surv"
    if (use_boot) {
      plotdata <- adj$boot_adjsurv
    } else {
      plotdata <- adj$adjsurv
    }
  } else if (inherits(adj, "adjustedcif")) {
    est_type <- "cif"
    if (use_boot) {
      plotdata <- adj$boot_adjcif
    } else {
      plotdata <- adj$adjcif
    }
  }

  # set groups
  if (is.null(group_1) | is.null(group_2)) {
    group_1 <- levels(plotdata$group)[1]
    group_2 <- levels(plotdata$group)[2]
  }

  # keep only data with the groups of interest
  plotdata <- plotdata[plotdata$group==group_1 | plotdata$group==group_2,]
  plotdata$group <- factor(plotdata$group, levels=c(group_1, group_2))

  # get difference + confidence intervals
  time_points <- sort(unique(plotdata$time))
  diff_curve <- difference_function(adj=plotdata, times=time_points,
                                    est=est_type, interpolation=interpolation,
                                    conf_int=conf_int, conf_level=conf_level)
  if (conf_int) {
    colnames(diff_curve) <- c("time", "diff", "se", "ci_lower", "ci_upper")
  } else {
    colnames(diff_curve) <- c("time", "diff")
  }

  # keep only some times
  if (!is.null(times)) {
    diff_curve$group <- factor("1")
    diff_curve <- specific_times(plotdata=diff_curve, times=times, est="diff",
                                 interpolation=interpolation)
    diff_curve$group <- NULL
  }

  # calculate p-values using a one-sample t-test
  if (conf_int) {
    t_stat <- (diff_curve$diff - 0) / diff_curve$se
    diff_curve$p_value <- 2 * stats::pnorm(abs(t_stat), lower.tail=FALSE)
  }

  return(diff_curve)
}
