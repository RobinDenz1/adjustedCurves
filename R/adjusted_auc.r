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
area_under_curve <- function(adj, from, to, conf_int, conf_level,
                             interpolation, subdivisions) {

  # silence checks
  . <- se <- group <- NULL

  mode <- ifelse(inherits(adj, "adjustedsurv"), "surv", "cif")
  boot_str <- paste0("boot_adj", mode)

  ## Using multiple imputation
  if (!is.null(adj$mids_analyses)) {

    # perform analysis on each MI object
    mids_out <- vector(mode="list", length=length(adj$mids_analyses))
    for (i in seq_len(length(adj$mids_analyses))) {

      mids_out[[i]] <- area_under_curve(adj$mids_analyses[[i]],
                                        to=to, from=from, conf_int=conf_int,
                                        conf_level=conf_level,
                                        interpolation=interpolation,
                                        subdivisions=subdivisions)
    }
    mids_out <- dplyr::bind_rows(mids_out)

    # pool results
    if (!conf_int) {
      out <- mids_out %>%
        dplyr::group_by(., group) %>%
        dplyr::summarise(auc=mean(auc))
      out <- as.data.frame(out)
    } else {
      out <- mids_out %>%
        dplyr::group_by(., group) %>%
        dplyr::summarise(auc=mean(auc),
                         se=mean(se),
                         n_boot=min(n_boot))

      auc_ci <- confint_surv(surv=out$auc,
                             se=out$se,
                             conf_level=conf_level,
                             conf_type="plain")
      out <- data.frame(group=out$group,
                        auc=out$auc,
                        se=out$se,
                        ci_lower=auc_ci$left,
                        ci_upper=auc_ci$right,
                        n_boot=out$n_boot)
    }

    return(out)

  ## single analysis
  } else {

    if (conf_int & !is.null(adj[boot_str])) {

      n_boot <- max(adj$boot_data$boot)
      booted_aucs <- vector(mode="list", length=n_boot)

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
        adj_auc <- area_under_curve(fake_object, from=from, to=to,
                                    conf_int=FALSE,
                                    interpolation=interpolation,
                                    subdivisions=subdivisions)
        adj_auc$boot <- i
        booted_aucs[[i]] <- adj_auc
      }
      booted_aucs <- as.data.frame(dplyr::bind_rows(booted_aucs))
    }

    if (inherits(adj, "adjustedsurv")) {
      levs <- levels(adj$adjsurv$group)
    } else {
      levs <- levels(adj$adjcif$group)
    }

    auc_vec <- vector(mode="numeric", length=length(levs))
    for (i in seq_len(length(levs))) {

      if (mode=="surv") {
        surv_dat <- adj$adjsurv[adj$adjsurv$group==levs[i], ]
      } else {
        surv_dat <- adj$adjcif[adj$adjcif$group==levs[i], ]
      }
      surv_dat$group <- NULL
      surv_dat$se <- NULL
      surv_dat$ci_lower <- NULL
      surv_dat$ci_upper <- NULL
      surv_dat$boot <- NULL

      if (interpolation=="steps") {
        auc <- exact_stepfun_integral(surv_dat, from=from, to=to, est=mode)
      } else if (interpolation=="linear") {
        times <- seq(from, to, (to-from)/subdivisions)
        diff_dat <- data.frame(times=times)
        diff_dat[, mode] <- vapply(X=times, FUN=read_from_linear_function,
                                   FUN.VALUE=numeric(1), data=surv_dat,
                                   est=mode)
        auc <- trapezoid_integral(x=diff_dat$time, y=diff_dat[, mode])
      }

      auc_vec[i] <- auc
    }

    out <- data.frame(group=levs, auc=auc_vec)

    if (conf_int & !is.null(adj[boot_str])) {

      out_boot <- booted_aucs %>%
        dplyr::group_by(., group) %>%
        dplyr::summarise(se=stats::sd(auc, na.rm=TRUE),
                         n_boot=sum(!is.na(auc)))
      out_boot$group <- NULL

      auc_ci <- confint_surv(surv=out$auc,
                             se=out_boot$se,
                             conf_level=conf_level,
                             conf_type="plain")
      out <- data.frame(group=out$group,
                        auc=out$auc,
                        se=out_boot$se,
                        ci_lower=auc_ci$left,
                        ci_upper=auc_ci$right,
                        n_boot=out_boot$n_boot)

    }
    return(out)
  }
}

## function to calculate the restricted mean survival time of each
## adjusted survival curve previously estimated using the adjustedsurv function
#' @export
adjusted_rmst <- function(adjsurv, to, from=0, conf_int=FALSE,
                          conf_level=0.95, interpolation="steps",
                          subdivisions=1000) {

  check_inputs_adj_rmst(adjsurv=adjsurv, from=from, to=to, conf_int=conf_int)

  # set to FALSE if it can't be done
  if (conf_int & is.null(adjsurv$boot_adjsurv)) {
    conf_int <- FALSE
  }

  out <- area_under_curve(adj=adjsurv, to=to, from=from,
                          conf_int=conf_int, conf_level=conf_level,
                          interpolation=interpolation,
                          subdivisions=subdivisions)

  if (conf_int) {
    colnames(out) <- c("group", "rmst", "se", "ci_lower", "ci_upper", "n_boot")
  } else {
    colnames(out) <- c("group", "rmst")
  }

  return(out)
}

## function to calculate the restricted mean time lost of each
## adjusted CIF previously estimated using the adjustedsurv/adjustedcif function
#' @export
adjusted_rmtl <- function(adj, to, from=0, conf_int=FALSE,
                          conf_level=0.95, interpolation="steps",
                          subdivisions=1000) {

  check_inputs_adj_rmtl(adj=adj, from=from, to=to, conf_int=conf_int)

  # set to FALSE if it can't be done
  if (conf_int & is.null(adj$boot_adjsurv) & is.null(adj$boot_adjcif)) {
    conf_int <- FALSE
  }

  # calculate area under curve
  out <- area_under_curve(adj=adj, to=to, from=from, conf_int=conf_int,
                          conf_level=conf_level,
                          interpolation=interpolation,
                          subdivisions=subdivisions)

  # take area above curve instead
  if (inherits(adj, "adjustedsurv")) {

    full_area <- (to - from) * 1
    out$auc <- full_area - out$auc

    # recalculate bootstrap CI
    if (conf_int) {
      auc_ci <- confint_surv(surv=out$auc,
                             se=out$se,
                             conf_level=conf_level,
                             conf_type="plain")
      out$ci_lower <- auc_ci$left
      out$ci_upper <- auc_ci$right
    }
  }

  if (conf_int) {
    colnames(out) <- c("group", "rmtl", "se", "ci_lower", "ci_upper", "n_boot")
  } else {
    colnames(out) <- c("group", "rmtl")
  }

  return(out)
}
