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
#' @export
adjusted_median_survival <- function(adjsurv, verbose=TRUE) {

  plotdata <- adjsurv$adjsurv
  plotdata <- plotdata[!is.na(plotdata$surv),]
  plotdata <- data.frame(surv_time=plotdata$time,
                         time=plotdata$surv,
                         group=plotdata$group)

  levs <- levels(plotdata$group)

  out <- vector(mode="numeric", length=length(levs))
  for (i in seq_len(length(levs))) {
    temp_dat <- plotdata[plotdata$group==levs[i],]
    temp_dat$group <- NULL

    out[[i]] <- read_from_step_function(0.5, step_data=temp_dat,
                                        est="surv_time")
  }

  out <- data.frame(group=levs, median_surv=out)

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
