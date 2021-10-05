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

## function to calculate the restricted mean survival time of each
## adjusted survival curve previously estimated using the adjustedsurv function
#' @export
adjusted_rmst <- function(adjsurv, to, from=0, use_boot=FALSE,
                          conf_level=0.95) {

  ## Using multiple imputation
  if (!is.null(adjsurv$mids_analyses)) {

    len <- length(adjsurv$mids_analyses)
    mids_out <- rmst_ests <- rmst_se <- n_boots <- vector(mode="list",
                                                          length=len)
    for (i in seq_len(len)) {

      results_imp <- adjusted_rmst(adjsurv$mids_analyses[[i]],
                                   to=to, from=from, use_boot=use_boot,
                                   conf_level=conf_level)
      mids_out[[i]] <- results_imp
      rmst_ests[[i]] <- results_imp$rmsts
      rmst_se[[i]] <- results_imp$rmsts_se
      n_boots[[i]] <- results_imp$n_boot

    }
    rmst_dat <- dplyr::bind_rows(rmst_ests)
    rmst_dat <- apply(rmst_dat, 2, mean)

    out <- list(mids_analyses=mids_out,
                rmsts=rmst_dat,
                from=from,
                to=to)
    class(out) <- "adjusted_rmst"

    if (use_boot) {

      rmst_se_dat <- dplyr::bind_rows(rmst_se)
      rmst_se_dat <- apply(rmst_se_dat, 2, mean)

      n_boots <- dplyr::bind_rows(n_boots)
      n_boots <- apply(n_boots, 2, mean)

      surv_ci <- confint_surv(surv=rmst_dat,
                              se=rmst_se_dat,
                              conf_level=conf_level,
                              conf_type="plain")

      out$conf_level <- conf_level
      out$n_boot <- n_boots
      out$rmsts_se <- rmst_se_dat
      out$rmsts_ci_lower <- surv_ci$left
      out$rmsts_ci_upper <- surv_ci$right

    }

    return(out)

  ## single analysis
  } else {

    check_inputs_adj_rmst(adjsurv=adjsurv, from=from, to=to, use_boot=use_boot)

    if (use_boot & !is.null(adjsurv$boot_adjsurv)) {

      n_boot <- max(adjsurv$boot_data$boot)
      booted_rmsts <- vector(mode="list", length=n_boot)

      for (i in 1:n_boot) {

        # select one bootstrap data set each
        boot_dat <- adjsurv$boot_data[adjsurv$boot_data$boot==i,]

        # create fake adjustedsurv object
        fake_adjsurv <- list(adjsurv=boot_dat)
        class(fake_adjsurv) <- "adjustedsurv"

        # recursion call
        adj_rmst <- adjusted_rmst(fake_adjsurv, from=from, to=to,
                                  use_boot=FALSE)

        booted_rmsts[[i]] <- adj_rmst$rmsts
      }
      booted_rmsts <- as.data.frame(dplyr::bind_rows(booted_rmsts))
    }

    levs <- unique(adjsurv$adjsurv$group)

    rmsts <- vector(mode="numeric", length=length(levs))
    for (i in seq_len(length(levs))) {

      surv_dat <- adjsurv$adjsurv[adjsurv$adjsurv$group==levs[i],]
      surv_dat$group <- NULL
      surv_dat$sd <- NULL
      surv_dat$se <- NULL
      surv_dat$ci_lower <- NULL
      surv_dat$ci_upper <- NULL
      surv_dat$boot <- NULL

      rmst <- exact_stepfun_integral(surv_dat, from=from, to=to)
      rmsts[i] <- rmst

    }
    names(rmsts) <- levs

    out <- list(rmsts=rmsts,
                from=from,
                to=to)
    class(out) <- "adjusted_rmst"

    if (use_boot & !is.null(adjsurv$boot_adjsurv)) {

      n_boot_rmst <- apply(booted_rmsts, 2, function(x){sum(!is.na(x))})
      names(n_boot_rmst) <- levs

      out$conf_level <- conf_level
      out$n_boot <- n_boot_rmst
      out$booted_rmsts <- booted_rmsts
      out$rmsts_se <- apply(booted_rmsts, 2, stats::sd, na.rm=TRUE)
      out$rmsts_ci_lower <- apply(booted_rmsts, 2, stats::quantile,
                                  probs=1-conf_level, na.rm=TRUE)
      out$rmsts_ci_upper <- apply(booted_rmsts, 2, stats::quantile,
                                  probs=conf_level, na.rm=TRUE)
    }

    return(out)

  }

}

## print method for the adjusted_rmst function
#' @export
print.adjusted_rmst <- function(x, digits=5, ...) {

  cat("------------------------------------------------------------------\n")
  cat("       Confounder-Adjusted Restricted Mean Survival Time\n")
  cat("------------------------------------------------------------------\n")
  cat("\n")
  cat("Using the interval:", x$from, "to", x$to, "\n")
  cat("\n")

  all_data <- rmst_data_frame(x)

  all_data <- round(all_data, digits)
  print(all_data)

  cat("------------------------------------------------------------------\n")

  # also silently return that data.frame
  return(invisible(all_data))
}

## plot method for the adjusted_rmst function
#' @importFrom rlang .data
#' @export
plot.adjusted_rmst <- function(x, draw_ci=TRUE, color=TRUE, point_size=2,
                               ci_size=1, ci_width=0.7, xlab="",
                               ylab="Adjusted RMST", title=NULL,
                               gg_theme=ggplot2::theme_classic(), ...) {

  # get data from object
  plotdata <- rmst_data_frame(x)
  colnames(plotdata) <- c("RMST", "RMST_SE", "ci_lower",
                          "ci_upper", "n_boot")
  plotdata$group <- rownames(plotdata)

  # plot it
  mapping <- ggplot2::aes(x=.data$group, y=.data$RMST, color=.data$group)

  if (!color) {
    mapping$colour <- NULL
  }

  if (!is.null(x$rmsts_ci_lower) & draw_ci) {

    p <- ggplot2::ggplot(plotdata, mapping) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=.data$ci_lower,
                                          ymax=.data$ci_upper),
                             size=ci_size, width=ci_width) +
      ggplot2::geom_point(size=point_size) +
      gg_theme +
      ggplot2::labs(x=xlab, y=ylab) +
      ggplot2::theme(legend.position="none")

  } else {

    p <- ggplot2::ggplot(plotdata, mapping) +
      ggplot2::geom_point(size=point_size) +
      gg_theme +
      ggplot2::labs(x=xlab, y=ylab) +
      ggplot2::theme(legend.position="none")

  }

  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }

  return(p)
}

## structure output of adjusted_rmst objects into neat data.frame
rmst_data_frame <- function(x) {

  if (!is.null(x$rmsts_ci_lower)) {
    all_data <- rbind(x$rmsts, x$rmsts_se, x$rmsts_ci_lower,
                      x$rmsts_ci_upper, x$n_boot)

    ci_lower_name <- paste0(round(x$conf_level*100, 2), "% CI (lower)")
    ci_upper_name <- paste0(round(x$conf_level*100, 2), "% CI (upper)")

    rownames(all_data) <- c("RMST", "RMST SE", ci_lower_name, ci_upper_name,
                            "N Boot")
    colnames(all_data) <- colnames(all_data)
    all_data <- t(all_data)
  } else {
    all_data <- data.frame(x$rmsts)
    colnames(all_data) <- c("RMST")
    rownames(all_data) <- rownames(all_data)
  }

  return(as.data.frame(all_data))
}
