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
adjusted_median_survival <- function(adjsurv, use_boot=FALSE,
                                     verbose=TRUE) {

  # define those to remove Notes in devtools::check()
  . <- time <- group <- surv <- NULL

  if (use_boot & is.null(adjsurv$boot_adjsurv)) {
    warning("Cannot use bootstrapped estimates as they were not estimated.",
            " Need bootstrap=TRUE in adjustedsurv() call.")
    plotdata <- adjsurv$adjsurv
  } else if (use_boot) {
    plotdata <- adjsurv$boot_adjsurv
  } else {
    plotdata <- adjsurv$adjsurv
  }

  out <- plotdata %>%
    dplyr::group_by(., group) %>%
    dplyr::summarise(median_surv=max(time[surv >= 0.5], na.rm=TRUE),
                     valid=ifelse(min(surv) <= 0.5, TRUE, FALSE),
                     .groups="drop_last")
  out$median_surv[!out$valid] <- NA
  out$valid <- NULL
  out <- as.data.frame(out)

  out_print <- out
  colnames(out_print) <- NULL

  if (verbose) {
    cat("----------------------------------\n")
    cat("Adjusted Median Survival Time\n")
    cat("----------------------------------\n")
    print(out_print, row.names=FALSE)
    cat("----------------------------------\n")
  }

  # also silently return that data.frame
  return(invisible(out))

}
