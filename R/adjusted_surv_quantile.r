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

## calculate confounder adjusted survival time quantile
#' @export
adjusted_surv_quantile <- function(adjsurv, p=0.5, conf_int=FALSE,
                                   use_boot=FALSE, interpolation="steps") {

  check_inputs_surv_q(adjsurv=adjsurv, conf_int=conf_int, p=p,
                      use_boot=use_boot, interpolation=interpolation)

  if (use_boot) {
    plotdata <- adjsurv$boot_adjsurv
  } else {
    plotdata <- adjsurv$adjsurv
  }

  # remove NAs
  plotdata <- plotdata[!is.na(plotdata$surv),]

  levs <- levels(plotdata$group)

  out <- vector(mode="list", length=length(levs))
  for (i in seq_len(length(levs))) {
    temp_dat <- plotdata[plotdata$group==levs[i],]
    temp_dat$group <- NULL

    val <- vapply(X=p, FUN=read_from_fun, FUN.VALUE=numeric(1), data=temp_dat,
                  est="time", time="surv", interpolation=interpolation)
    out_i <- data.frame(p=p, group=levs[i], q_surv=val)

    if (conf_int) {
      temp_dat <- temp_dat[!is.na(temp_dat$ci_lower) &
                           !is.na(temp_dat$ci_upper) ,]
      out_i$ci_lower <- vapply(X=p, FUN=read_from_fun, FUN.VALUE=numeric(1),
                               data=temp_dat, est="time", time="ci_lower",
                               interpolation=interpolation)
      out_i$ci_upper <- vapply(X=p, FUN=read_from_fun, FUN.VALUE=numeric(1),
                               data=temp_dat, est="time", time="ci_upper",
                               interpolation=interpolation)
    }
    out[[i]] <- out_i
  }
  out <- dplyr::bind_rows(out)

  return(out)
}
