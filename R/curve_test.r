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
#' @export
test_curve_equality <- function(adjsurv, to, from=0, conf_level=0.95) {

  check_inputs_adj_test(adjsurv=adjsurv, from=from, to=to)

  est <- ifelse(class(adjsurv)=="adjustedsurv", "surv", "cif")

  # just two treatments, standard procedure
  if (!adjsurv$categorical) {
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

    # shift bootstrap distribution
    diff_under_H0 <- stats_vec - mean(stats_vec)
    p_value <- mean(abs(diff_under_H0) > abs(observed_diff_integral))

    # bootstrap confidence interval
    conf_int <- stats::quantile(stats_vec,
                                probs=c((1-conf_level)/2, 1-((1-conf_level)/2)),
                                na.rm=T)
    names(conf_int) <- c("ci_lower", "ci_upper")

    # put together output
    out <- list(diff_curves=diff_curves,
                diff_integrals=stats_vec,
                observed_diff_curve=observed_diff_curve,
                observed_diff_integral=observed_diff_integral,
                p_value=p_value,
                n_boot=length(stats_vec),
                from=from,
                to=to,
                kind=est,
                conf_int=conf_int,
                conf_level=conf_level,
                categorical=F,
                treat_labs=unique(adjsurv$adjsurv$group))

  ## more than one treatments -> perform pairwise comparisons
  } else {
    combs <- all_combs_length_2(unique(adjsurv$adjsurv$group))

    out <- list()
    for (i in 1:length(combs)) {

      # get first and second group
      group_0 <- substring(combs[i], 1, 1)
      group_1 <- substring(combs[i], 2, 2)

      # create pseudo adjustedsurv, adjustedcif object
      if (class(adjsurv)=="adjustedsurv") {

        observed_dat <- adjsurv$adjsurv[which(adjsurv$adjsurv$group %in%
                                                c(group_0, group_1)),]
        boot_dat <- adjsurv$boot_data[which(adjsurv$boot_data$group %in%
                                              c(group_0, group_1)),]
        fake_adjsurv <- list(adjsurv=observed_dat,
                             boot_data=boot_dat,
                             categorical=F)
        class(fake_adjsurv) <- "adjustedsurv"

      } else {

        observed_dat <- adjsurv$adjcif[which(adjsurv$adjcif$group %in%
                                               c(group_0, group_1)),]
        boot_dat <- adjsurv$boot_data[which(adjsurv$boot_data$group %in%
                                              c(group_0, group_1)),]
        fake_adjsurv <- list(adjcif=observed_dat,
                             boot_data=boot_dat,
                             categorical=F)
        class(fake_adjsurv) <- "adjustedcif"
      }

      # recursion call
      pair <- test_curve_equality(adjsurv=fake_adjsurv,
                                  from=from,
                                  to=to,
                                  conf_level=conf_level)
      out[[paste(group_0, group_1)]] <- pair
    }
    out$categorical <- T
  }
  class(out) <- "curve_test"

  return(out)
}

## print method for curve_test objects
#' @export
print.curve_test <- function(x, ...) {

  if (!x$categorical) {
    if(x$kind=="surv") {
      title <- "Pepe-Flemming Test of Equality of Two Adjusted Survival Curves \n"
    } else {
      title <- "Pepe-Flemming Test of Equality of Two Adjusted CIFs \n"
    }

    cat("------------------------------------------------------------------\n")
    cat(title)
    cat("------------------------------------------------------------------\n")
    cat("Group =", x$treat_labs[1], "vs. Group =", x$treat_labs[2], "\n")
    cat("The equality was tested for the time interval:", x$from, "to", x$to, "\n")
    cat("Observed Integral of the difference:", x$observed_diff_integral, "\n")
    cat("Bootstrap standard deviation:", stats::sd(x$diff_integrals), "\n")
    cat(x$conf_level*100, "% bootstrap confidence interval: [",
        x$conf_int["ci_lower"], ", ", x$conf_int["ci_upper"], "]\n", sep="")
    cat("P-Value:", x$p_value, "\n\n")
    cat("Calculated using", x$n_boot, "bootstrap replications.\n")
    cat("------------------------------------------------------------------\n")
  } else {
    for (i in 1:(length(x)-1))  {
      print.curve_test(x=x[[i]])
    }
  }
}

## function that calculates all possible combinations
## of length 2 without replacement and without order
all_combs_length_2 <- function(treat_labs) {

  combs <- list()
  for (i in 1:length(treat_labs)) {
    for (j in 1:length(treat_labs)) {

      # no cases with the same group twice and order does not matter
      if (i != j & !paste0(treat_labs[j], treat_labs[i]) %in% combs) {
        combs[[length(combs)+1]] <- paste0(treat_labs[i], treat_labs[j])
      }
    }
  }
  combs <- unlist(combs)
  return(combs)
}
