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

## calculate confounder adjusted median survival times
#' @importFrom dplyr %>%
#' @export
adjusted_median_survival <- function(adjsurv, use_boot=F, conf_level=0.95,
                                     verbose=T) {

  # define those to remove Notes in devtools::check()
  . <- i <- time <- group <- median_surv_i <- boot <- surv <- NULL

  if (use_boot & is.null(adjsurv$boot_data)) {

    stop("Cannot use bootstrapped estimates as they were not estimated.",
         " Need bootstrap=TRUE in adjustedsurv() call.")

  } else if (use_boot) {

    out <- adjsurv$boot_data %>%
      dplyr::group_by(., boot, group) %>%
      dplyr::summarise(median_surv_i=max(time[surv >= 0.5], na.rm=T),
                       .groups="drop_last") %>%
      dplyr::group_by(., group) %>%
      dplyr::summarise(median_surv=mean(median_surv_i, na.rm=T),
                       ci_lower=stats::quantile(median_surv_i,
                                                probs=(1-conf_level)/2,
                                                na.rm=T),
                       ci_upper=stats::quantile(median_surv_i,
                                                probs=1-((1-conf_level)/2),
                                                na.rm=T),
                       n=sum(!is.na(median_surv_i)),
                       .groups="drop_last")
    out <- as.data.frame(out)

  } else {

    out <- adjsurv$adjsurv %>%
      dplyr::group_by(., group) %>%
      dplyr::summarise(median_surv=max(time[surv >= 0.5], na.rm=T),
                       .groups="drop_last")
    out <- as.data.frame(out)

  }

  if (verbose) {
    cat("---------------------------------------------------\n")
    cat("Adjusted Median Survival Time\n")
    cat("---------------------------------------------------\n")
    print(out, row.names=F)
    cat("---------------------------------------------------\n")
  }

  # also silently return that data.frame
  return(invisible(out))
}
