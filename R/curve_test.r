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

## Hypothesis-Test for the equality of two adjusted survival curves
# TODO:
# - should work with pairwise comparisons
# - give confidence interval for difference integral
#' @export
test_curve_equality <- function(adjsurv, to, from=0) {

  check_inputs_adj_test(adjsurv=adjsurv, from=from, to=to)

  est <- ifelse(class(adjsurv)=="adjustedsurv", "surv", "cif")

  # calculate the integral of the difference for every bootstrap sample
  stats_vec <- vector(mode="numeric", length=max(adjsurv$boot_data$boot))
  curve_list <- vector(mode="list", length=max(adjsurv$boot_data$boot))

  for (i in 1:max(adjsurv$boot_data$boot)) {

    # select one bootstrap data set each
    boot_dat <- adjsurv$boot_data[adjsurv$boot_data$boot==i,]

    # every relevant point in time
    times <- sort(unique(boot_dat$time))

    # new curve of the difference
    surv_diff <- exact_stepfun_difference(adjsurv=boot_dat, times=times,
                                          est=est)
    curve_list[[i]] <- surv_diff

    # integral of that curve
    diff_integral <- exact_stepfun_integral(surv_diff, to=to, from=from,
                                            est=est)
    stats_vec[i] <- diff_integral
  }

  diff_curves <- as.data.frame(dplyr::bind_rows(curve_list))

  ## use those bootstrapped integrals for the calculation of
  ## a p-value, by shifting the observed distribution to 0 and
  ## comparing it to the actually observed value

  # actually observed values
  if (est=="surv") {
    times <- sort(unique(adjsurv$adjsurv$time))
    observed_diff_curve <- exact_stepfun_difference(adjsurv=adjsurv$adjsurv,
                                                    times=times, est=est)
  } else {
    times <- sort(unique(adjsurv$adjcif$time))
    observed_diff_curve <- exact_stepfun_difference(adjsurv=adjsurv$adjcif,
                                                    times=times, est=est)
  }

  observed_diff_integral <- exact_stepfun_integral(observed_diff_curve,
                                                   to=to, from=from,
                                                   est=est)

  # remove NA values
  stats_vec <- stats_vec[!is.na(stats_vec)]

  # shit bootstrap distribution
  diff_under_H0 <- stats_vec - mean(stats_vec)
  p_value <- mean(abs(diff_under_H0) > abs(observed_diff_integral))

  # put together output
  out <- list(diff_curves=diff_curves,
              diff_integrals=stats_vec,
              observed_diff_curve=observed_diff_curve,
              observed_diff_integral=observed_diff_integral,
              p_value=p_value,
              n_boot=length(stats_vec),
              from=from,
              to=to,
              kind=est)
  class(out) <- "curve_test"

  return(out)
}

## print method for curve_test objects
#' @export
print.curve_test <- function(x, ...) {

  if(x$kind=="surv") {
    title <- "Pepe-Flemming Test of Equality of Two Adjusted Survival Curves \n"
  } else {
    title <- "Pepe-Flemming Test of Equality of Two Adjusted CIFs \n"
  }

  cat("------------------------------------------------------------------\n")
  cat(title)
  cat("------------------------------------------------------------------\n")
  cat("The equality was tested for the time interval:", x$from, "to", x$to, "\n")
  cat("Observed Integral of the difference:", x$observed_diff_integral, "\n")
  cat("Bootstrap standard deviation:", stats::sd(x$diff_integrals), "\n")
  cat("P-Value:", x$p_value, "\n\n")
  cat("Calculated using", x$n_boot, "bootstrap replications.\n")
  cat("------------------------------------------------------------------\n")
}
