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

## calculate area under curve of adjustedsurv, adjustedcif objects
area_under_curve <- function(adj, from, to, use_boot, conf_level) {

  mode <- ifelse(inherits(adj, "adjustedsurv"), "surv", "cif")
  boot_str <- paste0("boot_adj", mode)

  ## Using multiple imputation
  if (!is.null(adj$mids_analyses)) {

    len <- length(adj$mids_analyses)
    mids_out <- rmst_ests <- rmst_se <- n_boots <- vector(mode="list",
                                                          length=len)
    for (i in seq_len(len)) {

      results_imp <- area_under_curve(adj$mids_analyses[[i]],
                                      to=to, from=from, use_boot=use_boot,
                                      conf_level=conf_level)
      mids_out[[i]] <- results_imp
      rmst_ests[[i]] <- results_imp$auc
      rmst_se[[i]] <- results_imp$auc_se
      n_boots[[i]] <- results_imp$n_boot

    }
    rmst_dat <- dplyr::bind_rows(rmst_ests)
    rmst_dat <- apply(rmst_dat, 2, mean)

    out <- list(mids_analyses=mids_out,
                auc=rmst_dat,
                from=from,
                to=to)

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
      out$auc_se <- rmst_se_dat
      out$auc_ci_lower <- surv_ci$left
      out$auc_ci_upper <- surv_ci$right

    }

    return(out)

  ## single analysis
  } else {

    if (use_boot & !is.null(adj[boot_str])) {

      n_boot <- max(adj$boot_data$boot)
      booted_rmsts <- vector(mode="list", length=n_boot)

      for (i in seq_len(n_boot)) {

        # select one bootstrap data set each
        boot_dat <- adj$boot_data[adj$boot_data$boot==i, ]

        # create fake adjustedsurv object
        if (mode=="surv") {
          fake_object <- list(adjsurv=boot_dat)
          class(fake_object) <- "adjustedsurv"
        } else {
          fake_object <- list(adjcif=boot_dat)
          class(fake_object) <- "adjustedcif"
        }

        # recursion call
        adj_rmst <- area_under_curve(fake_object, from=from, to=to,
                                     use_boot=FALSE)

        booted_rmsts[[i]] <- adj_rmst$auc
      }
      booted_rmsts <- as.data.frame(dplyr::bind_rows(booted_rmsts))
    }

    if (inherits(adj, "adjustedsurv")) {
      levs <- unique(adj$adjsurv$group)
    } else {
      levs <- unique(adj$adjcif$group)
    }

    rmsts <- vector(mode="numeric", length=length(levs))
    for (i in seq_len(length(levs))) {

      if (mode=="surv") {
        surv_dat <- adj$adjsurv[adj$adjsurv$group==levs[i], ]
      } else {
        surv_dat <- adj$adjcif[adj$adjcif$group==levs[i], ]
      }
      surv_dat$group <- NULL
      surv_dat$sd <- NULL
      surv_dat$se <- NULL
      surv_dat$ci_lower <- NULL
      surv_dat$ci_upper <- NULL
      surv_dat$boot <- NULL

      rmst <- exact_stepfun_integral(surv_dat, from=from, to=to, est=mode)
      rmsts[i] <- rmst

    }
    names(rmsts) <- levs

    out <- list(auc=rmsts,
                from=from,
                to=to)

    if (use_boot & !is.null(adj[boot_str])) {

      n_boot_rmst <- apply(booted_rmsts, 2, function(x) {sum(!is.na(x))})
      names(n_boot_rmst) <- levs

      out$conf_level <- conf_level
      out$n_boot <- n_boot_rmst
      out$booted_auc <- booted_rmsts
      out$auc_se <- apply(booted_rmsts, 2, stats::sd, na.rm=TRUE)
      out$auc_ci_lower <- apply(booted_rmsts, 2, stats::quantile,
                                probs=1-conf_level, na.rm=TRUE)
      out$auc_ci_upper <- apply(booted_rmsts, 2, stats::quantile,
                                probs=conf_level, na.rm=TRUE)
    }

    return(out)
  }
}

## function to calculate the restricted mean survival time of each
## adjusted survival curve previously estimated using the adjustedsurv function
#' @export
adjusted_rmst <- function(adjsurv, to, from=0, use_boot=FALSE,
                          conf_level=0.95) {

  check_inputs_adj_rmst(adjsurv=adjsurv, from=from, to=to, use_boot=use_boot)

  # set to FALSE if it can't be done
  if (use_boot & is.null(adjsurv$boot_adjsurv)) {
    use_boot <- FALSE
  }

  out <- area_under_curve(adj=adjsurv, to=to, from=from,
                          use_boot=use_boot, conf_level=conf_level)
  class(out) <- "adjusted_rmst"
  return(out)
}

## function to calculate the restricted mean time lost of each
## adjusted CIF previously estimated using the adjustedsurv/adjustedcif function
#' @export
adjusted_rmtl <- function(adj, to, from=0, use_boot=FALSE,
                          conf_level=0.95) {

  check_inputs_adj_rmtl(adj=adj, from=from, to=to, use_boot=use_boot)

  # set to FALSE if it can't be done
  if (use_boot & is.null(adj$boot_adjsurv) & is.null(adj$boot_adjcif)) {
    use_boot <- FALSE
  }

  # calculate area under curve
  out <- area_under_curve(adj=adj, to=to, from=from, use_boot=use_boot,
                          conf_level=conf_level)
  class(out) <- "adjusted_rmtl"

  # take area above curve instead
  if (inherits(adj, "adjustedsurv")) {

    full_area <- (to - from) * 1
    out$auc <- full_area - out$auc

    # recalculate SE / CI
    if (use_boot) {
      out$booted_auc <- full_area - out$booted_auc
      out$auc_se <- apply(out$booted_auc, 2, stats::sd, na.rm=TRUE)
      out$auc_ci_lower <- apply(out$booted_auc, 2, stats::quantile,
                                probs=1-conf_level, na.rm=TRUE)
      out$auc_ci_upper <- apply(out$booted_auc, 2, stats::quantile,
                                probs=conf_level, na.rm=TRUE)
    }
  }
  return(out)
}

## print method for the adjusted_rmst function
#' @export
print.adjusted_rmst <- function(x, digits=5, ...) {

  cat("------------------------------------------------------------------\n")
  if (inherits(x, "adjusted_rmst")) {
    cat("       Confounder-Adjusted Restricted Mean Survival Time\n")
  } else {
    cat("       Confounder-Adjusted Restricted Mean Time Lost\n")
  }
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

## print method for the adjusted_rmtl function
#' @export
print.adjusted_rmtl <- function(x, digits=5, ...) {
  print.adjusted_rmst(x=x, digits=digits, ...)
}

## summary method for the adjusted_rmst function
#' @export
summary.adjusted_rmst <- function(object, digits=5, ...) {
  print.adjusted_rmst(x=object, digits=digits, ...)
}

## summary method for the adjusted_rmtl function
#' @export
summary.adjusted_rmtl <- function(object, digits=5, ...) {
  print.adjusted_rmst(x=object, digits=digits, ...)
}

## plot method for the adjusted_rmst function
#' @importFrom rlang .data
#' @export
plot.adjusted_rmst <- function(x, conf_int=TRUE, color=TRUE, point_size=2,
                               ci_size=1, ci_width=0.7, xlab="",
                               ylab="Adjusted RMST", title=NULL,
                               gg_theme=ggplot2::theme_classic(), ...) {
  requireNamespace("ggplot2")

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

  if (!is.null(x$auc_ci_lower) & conf_int) {

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

## plot method for the adjusted_rmst function
#' @importFrom rlang .data
#' @export
plot.adjusted_rmtl <- function(x, conf_int=TRUE, color=TRUE, point_size=2,
                               ci_size=1, ci_width=0.7, xlab="",
                               ylab="Adjusted RMTL", title=NULL,
                               gg_theme=ggplot2::theme_classic(), ...) {
  plot.adjusted_rmst(x=x, conf_int=conf_int, color=color,
                     point_size=point_size, ci_size=ci_size, ci_width=ci_width,
                     xlab=xlab, ylab=ylab, title=title, gg_theme=gg_theme,
                     ...)
}

## structure output of adjusted_rmst objects into neat data.frame
rmst_data_frame <- function(x) {

  if (!is.null(x$auc_ci_lower)) {
    all_data <- rbind(x$auc, x$auc_se, x$auc_ci_lower,
                      x$auc_ci_upper, x$n_boot)

    ci_lower_name <- paste0(round(x$conf_level*100, 2), "% CI (lower)")
    ci_upper_name <- paste0(round(x$conf_level*100, 2), "% CI (upper)")

    if (inherits(x, "adjusted_rmst")) {
      rownames(all_data) <- c("RMST", "RMST SE", ci_lower_name, ci_upper_name,
                              "N Boot")
    } else {
      rownames(all_data) <- c("RMTL", "RMTL SE", ci_lower_name, ci_upper_name,
                              "N Boot")
    }
    colnames(all_data) <- colnames(all_data)
    all_data <- t(all_data)
  } else {
    all_data <- data.frame(x$auc)
    if (inherits(x, "adjusted_rmst")) {
      colnames(all_data) <- c("RMST")
    } else {
      colnames(all_data) <- c("RMTL")
    }
    rownames(all_data) <- rownames(all_data)
  }

  return(as.data.frame(all_data))
}
