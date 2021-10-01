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
## or two adjusted cumulative incidence functions
#' @export
test_curve_equality <- function(adjsurv, to, from=0, conf_level=0.95) {

  # silence devtools::check() notes
  . <- comparison <- area_est <- p_val <- n_boot <- NULL

  est <- ifelse(class(adjsurv)=="adjustedsurv", "surv", "cif")

  if (est=="surv") {
    treat_labs <- unique(adjsurv$adjsurv$group)
  } else {
    treat_labs <- unique(adjsurv$adjcif$group)
  }

  ## using multiply imputed results
  if (!is.null(adjsurv$mids_analyses)) {

    if (adjsurv$categorical) {

      # call function once on every adjsurv object, extract values
      len <- length(adjsurv$mids_analyses)
      mids_out <- dat <- vector(mode="list", length=len)
      for (i in 1:len) {

        results_imp <- test_curve_equality(adjsurv$mids_analyses[[i]],
                                           to=to, from=from,
                                           conf_level=conf_level)
        mids_out[[i]] <- results_imp
        comp_names <- names(results_imp)

        # for each pairwise comparison, extract values
        for (j in 1:(length(results_imp)-1)) {

          row <- data.frame(area_est=results_imp[[j]]$observed_diff_integral,
                            area_se=results_imp[[j]]$integral_se,
                            p_val=results_imp[[j]]$p_value,
                            n_boot=results_imp[[j]]$n_boot,
                            comparison=comp_names[j])
          dat[[length(dat)+1]] <- row
        }

      }
      dat <- dplyr::bind_rows(dat)

      # pool values
      pooled_dat <- dat %>%
        dplyr::group_by(., comparison) %>%
        dplyr::summarise(area_est=mean(area_est),
                         area_se=mean(area_se),
                         p_value=pool_p_values(p_val),
                         n_boot=min(n_boot))

      # create one curve_test object for each comparison
      output <- list(mids_analyses=mids_out)
      for (i in 1:nrow(pooled_dat)) {

        pooled_ci <- confint_surv(surv=dat$area_est[i],
                                  se=dat$area_se[i],
                                  conf_level=conf_level,
                                  conf_type="plain")
        pooled_ci <- c(pooled_ci$left, pooled_ci$right)
        names(pooled_ci) <- c("ci_lower", "ci_upper")

        # do this so it shows up in print function
        fun_call <- match.call()
        fun_call$from <- from
        fun_call$to <- to

        out <- list(mids_p_values=dat$p_vals,
                    observed_diff_integral=dat$area_est[i],
                    p_value=dat$p_val[i],
                    n_boot=dat$n_boot[i],
                    kind=est,
                    integral_se=dat$area_se[i],
                    conf_int=pooled_ci,
                    categorical=F,
                    treat_labs=mids_out[[1]][[i]]$treat_labs,
                    call=fun_call)
        class(out) <- "curve_test"
        output[[dat$comparison[i]]] <- out

      }
      output$categorical <- T
      class(output) <- "curve_test"

      return(output)

    } else {

      len <- length(adjsurv$mids_analyses)
      mids_out <- vector(mode="list", length=len)
      area_ests <- area_se <- p_vals <- n_boots <- vector(mode="numeric",
                                                          length=len)
      for (i in 1:len) {

        results_imp <- test_curve_equality(adjsurv$mids_analyses[[i]],
                                           to=to, from=from,
                                           conf_level=conf_level)
        mids_out[[i]] <- results_imp
        area_ests[i] <- results_imp$observed_diff_integral
        area_se[i] <- results_imp$integral_se
        p_vals[i] <- results_imp$p_value
        n_boots[i] <- results_imp$n_boot

      }

      # pool the values
      observed_diff_integral <- mean(area_ests)
      pooled_se <- mean(area_se)
      pooled_ci <- confint_surv(surv=observed_diff_integral,
                                se=pooled_se,
                                conf_level=conf_level,
                                conf_type="plain")
      pooled_ci <- c(pooled_ci$left, pooled_ci$right)
      names(pooled_ci) <- c("ci_lower", "ci_upper")

      pooled_p_value <- pool_p_values(p_values=p_vals)

      # do this so it shows up in print function
      fun_call <- match.call()
      fun_call$from <- from
      fun_call$to <- to

      out <- list(mids_analyses=mids_out,
                  mids_p_values=p_vals,
                  observed_diff_integral=observed_diff_integral,
                  p_value=pooled_p_value,
                  n_boot=min(n_boots),
                  kind=est,
                  integral_se=pooled_se,
                  conf_int=pooled_ci,
                  categorical=F,
                  treat_labs=treat_labs,
                  call=fun_call)
      class(out) <- "curve_test"

    }

    return(out)

  ## using regular results
  } else {

    check_inputs_adj_test(adjsurv=adjsurv, from=from, to=to)

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

        # integral of that curve
        diff_integral <- exact_stepfun_integral(surv_diff, to=to, from=from,
                                                est=est)
        stats_vec[i] <- diff_integral

        # also append difference curve to list
        surv_diff$boot <- i
        curve_list[[i]] <- surv_diff
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
                  integral_se=stats::sd(stats_vec),
                  p_value=p_value,
                  n_boot=length(stats_vec),
                  kind=est,
                  conf_int=conf_int,
                  categorical=F,
                  treat_labs=treat_labs,
                  call=match.call())

      ## more than one treatments -> perform pairwise comparisons
    } else {

      if (class(adjsurv)=="adjustedsurv") {
        combs <- all_combs_length_2(unique(adjsurv$adjsurv$group))
      } else {
        combs <- all_combs_length_2(unique(adjsurv$adjcif$group))
      }


      out <- list()
      for (i in 1:length(combs)) {

        # get first and second group
        group_0 <- strsplit(combs[i], "\t")[[1]][1]
        group_1 <- strsplit(combs[i], "\t")[[1]][2]

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

}

## print method for curve_test objects
#' @export
print.curve_test <- function(x, ...) {

  if (!is.null(x$mids_analyses) & x$categorical) {

    for (i in 2:(length(x)-1))  {
      print.curve_test(x=x[[i]])
    }

  } else if (!x$categorical) {
    if(x$kind=="surv") {
      title <- "Pepe-Flemming Test of Equality of Two Adjusted Survival Curves \n"
    } else {
      title <- "Pepe-Flemming Test of Equality of Two Adjusted CIFs \n"
    }

    call_conf <- x$call$conf_level
    if (is.null(call_conf) || call_conf=="conf_level") {
      conf_level <- 95
    } else {
      conf_level <- call_conf * 100
    }

    cat("------------------------------------------------------------------\n")
    cat(title)
    cat("------------------------------------------------------------------\n")
    cat("Group =", toString(x$treat_labs[1]), "vs. Group =",
        toString(x$treat_labs[2]), "\n")
    cat("The equality was tested for the time interval:", x$call$from, "to",
        x$call$to, "\n")
    cat("Observed Integral of the difference:", x$observed_diff_integral, "\n")
    cat("Bootstrap standard error:", x$integral_se, "\n")
    cat(conf_level, "% bootstrap confidence interval: [",
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

## plot method for curve_test objects
#' @export
plot.curve_test <- function(x, type="curves", xlab=NULL, ylab=NULL,
                            title=NULL, ...) {

  # to remove devtools::check() Notes
  time <- surv <- boot <- integral <- NULL

  if (!is.null(x$mids_analyses)) {

    stop("There is no plot method for 'curve_test' objects fitted using",
         " multiply imputed datasets.")

  } else if (x$categorical) {

    comps <- names(x)

    if (type=="integral") {

      observed_diff_integrals <- list()
      diff_integrals <- list()
      for (i in 1:(length(x)-1)) {

        # observed diff curves
        obs_diff <- x[[i]]$observed_diff_integral
        observed_diff_integrals[[i]] <- data.frame(integral=obs_diff,
                                                   comp=comps[i])

        # bootstrapped diff curves
        diff <- x[[i]]$diff_integrals

        # remove NA values
        diff <- diff[!is.na(diff)]

        # shift bootstrap distribution
        diff <- diff - mean(diff)

        diff_integrals[[i]] <- data.frame(integral=diff,
                                          comp=comps[i])
      }
      observed_diff_integrals <- dplyr::bind_rows(observed_diff_integrals)
      diff_integrals <- dplyr::bind_rows(diff_integrals)

      if (is.null(xlab)) {
        xlab <- "Integrals of the difference under H0"
      }

      if (is.null(ylab)) {
        ylab <- "Density"
      }

      # plot it
      p <- ggplot2::ggplot(diff_integrals, ggplot2::aes(x=integral)) +
        ggplot2::geom_density() +
        ggplot2::geom_vline(xintercept=0, linetype="dashed") +
        ggplot2::geom_vline(data=observed_diff_integrals,
                            ggplot2::aes(xintercept=integral),
                            color="red") +
        ggplot2::theme_bw() +
        ggplot2::labs(x=xlab, y=ylab) +
        ggplot2::facet_wrap(~comp, scales="free")

    } else if (type=="curves") {

      if (is.null(ylab)) {
        ylab <- ifelse(x[[1]]$kind=="surv", "Difference in Survival",
                       "Difference in Cumulative Incidence")
      }

      if (is.null(xlab)) {
        xlab <- "Time"
      }

      observed_diff_curves <- list()
      diff_curves <- list()
      for (i in 1:(length(x)-1)) {

        # observed diff curves
        obs_diff <- x[[i]]$observed_diff_curve
        obs_diff$comp <- comps[i]
        observed_diff_curves[[i]] <- obs_diff

        # bootstrapped diff curves
        diff <- x[[i]]$diff_curves
        diff$comp <- comps[i]
        diff_curves[[i]] <- diff
      }
      observed_diff_curves <- dplyr::bind_rows(observed_diff_curves)
      diff_curves <- dplyr::bind_rows(diff_curves)

      p <- ggplot2::ggplot(diff_curves, ggplot2::aes(x=time, y=surv)) +
        ggplot2::geom_step(ggplot2::aes(group=boot), color="grey", alpha=0.8) +
        ggplot2::geom_step(data=observed_diff_curves, ggplot2::aes(x=time, y=surv)) +
        ggplot2::geom_hline(yintercept=0, linetype="dashed") +
        ggplot2::theme_bw() +
        ggplot2::labs(x=xlab, y=ylab) +
        ggplot2::facet_wrap(~comp, scales="free")
    }

  } else {

    if (type=="integral") {

      # remove NA values
      stats_vec <- x$diff_integrals[!is.na(x$diff_integrals)]

      # shift bootstrap distribution
      diff_under_H0 <- stats_vec - mean(stats_vec)

      if (is.null(xlab)) {
        xlab <- "Integrals of the difference under H0"
      }

      if (is.null(ylab)) {
        ylab <- "Density"
      }

      p <- ggplot2::ggplot(NULL, ggplot2::aes(x=diff_under_H0)) +
        ggplot2::geom_density() +
        ggplot2::geom_vline(xintercept=0, linetype="dashed") +
        ggplot2::geom_vline(xintercept=x$observed_diff_integral, color="red") +
        ggplot2::theme_bw() +
        ggplot2::labs(x=xlab, y=ylab)

    } else if (type=="curves") {

      if (is.null(ylab)) {
        ylab <- ifelse(x$kind=="surv", "Difference in Survival",
                       "Difference in Cumulative Incidence")
      }

      if (is.null(xlab)) {
        xlab <- "Time"
      }

      p <- ggplot2::ggplot(x$diff_curves, ggplot2::aes(x=time, y=surv)) +
        ggplot2::geom_step(ggplot2::aes(group=boot), color="grey", alpha=0.8) +
        ggplot2::geom_step(data=x$observed_diff_curve,
                           ggplot2::aes(x=time, y=surv)) +
        ggplot2::geom_hline(yintercept=0, linetype="dashed") +
        ggplot2::theme_bw() +
        ggplot2::labs(x=xlab, y=ylab)

    }
  }

  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }

  return(p)
}

## function that calculates all possible combinations
## of length 2 without replacement and without order
all_combs_length_2 <- function(treat_labs) {

  combs <- list()
  for (i in 1:length(treat_labs)) {
    for (j in 1:length(treat_labs)) {

      # no cases with the same group twice and order does not matter
      if (i != j & !paste(treat_labs[j], treat_labs[i], sep="\t") %in% combs) {
        combs[[length(combs)+1]] <- paste(treat_labs[i], treat_labs[j], sep="\t")
      }
    }
  }
  combs <- unlist(combs)
  return(combs)
}

## directly pool the p-values
pool_p_values <- function(p_values, method="licht_rubin", tol=1e-14) {

  # if they are all equal, just return that value
  if (stats::var(p_values)==0) {
    return(p_values[1])
  }

  # apply tolerance to p values that are 0
  p_values[p_values==0] <- p_values[p_values==0] + tol

  if (method=="licht_rubin") {

    # transform to z-scale
    z <- stats::qnorm(p_values)
    num <- mean(z)
    den <- sqrt(1 + stats::var(z))
    pooled <- stats::pnorm(num / den)

  } else if (method=="median") {

    pooled <- stats::median(p_values)

  } else if (method=="mean") {

    pooled <- mean(p_values)

  }

  return(pooled)

}
