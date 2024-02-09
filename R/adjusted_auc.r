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
                             interpolation) {

  # silence checks
  . <- se <- group <- NULL

  if (inherits(adj, "adjustedsurv")) {
    levs <- levels(adj$adjsurv$group)
  } else {
    levs <- levels(adj$adjcif$group)
  }

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
                                        interpolation=interpolation)
    }
    mids_out <- dplyr::bind_rows(mids_out)

    # pool results
    if (!conf_int) {
      out <- mids_out %>%
        dplyr::group_by(., to, group) %>%
        dplyr::summarise(auc=mean(auc))
      out <- as.data.frame(out)
    } else {
      out <- mids_out %>%
        dplyr::group_by(., to, group) %>%
        dplyr::summarise(auc=mean(auc),
                         se=mean(se),
                         n_boot=min(n_boot))

      auc_ci <- confint_surv(surv=out$auc,
                             se=out$se,
                             conf_level=conf_level,
                             conf_type="plain")
      out <- data.frame(to=out$to,
                        group=factor(out$group, levels=levs),
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
                                    interpolation=interpolation)
        adj_auc$boot <- i
        booted_aucs[[i]] <- adj_auc
      }
      booted_aucs <- as.data.frame(dplyr::bind_rows(booted_aucs))
    }

    ## core of the function
    out <- vector(mode="list", length=length(levs))
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

      auc <- exact_integral(data=surv_dat, from=from, to=to, est=mode,
                            interpolation=interpolation)
      out[[i]] <- data.frame(to=to, group=levs[i], auc=auc)
    }
    out <- dplyr::bind_rows(out)
    out$group <- factor(out$group, levels=levs)
    out <- out[order(out$to, out$group),]

    if (conf_int & !is.null(adj[boot_str])) {

      out_boot <- booted_aucs %>%
        dplyr::group_by(., to, group) %>%
        dplyr::summarise(se=stats::sd(auc, na.rm=TRUE),
                         n_boot=sum(!is.na(auc)))
      out_boot$group <- NULL

      auc_ci <- confint_surv(surv=out$auc,
                             se=out_boot$se,
                             conf_level=conf_level,
                             conf_type="plain")
      out <- data.frame(to=out$to,
                        group=factor(out$group, levels=levs),
                        auc=out$auc,
                        se=out_boot$se,
                        ci_lower=auc_ci$left,
                        ci_upper=auc_ci$right,
                        n_boot=out_boot$n_boot)

    }
    return(out)
  }
}

## get difference of AUC values + confidence intervals and p_value
auc_difference <- function(data, group_1, group_2, conf_int, conf_level) {

  if (is.null(group_1) | is.null(group_2)) {
    group_1 <- levels(data$group)[1]
    group_2 <- levels(data$group)[2]
  }

  dat_1 <- data[data$group==group_1, ]
  dat_2 <- data[data$group==group_2, ]

  out <- data.frame(to=dat_1$to, diff=dat_1$auc - dat_2$auc)

  if (conf_int) {
    out$se <- sqrt(dat_1$se^2 + dat_2$se^2)
    diff_ci <- confint_surv(surv=out$diff, se=out$se, conf_level=conf_level,
                            conf_type="plain")
    out$ci_lower <- diff_ci$left
    out$ci_upper <- diff_ci$right
    t_stat <- (out$diff - 0) / out$se
    out$p_value <- 2 * stats::pnorm(abs(t_stat), lower.tail=FALSE)
  }
  return(out)
}

## get ratio of AUC values + confidence intervals and p_value
auc_ratio <- function(data, group_1, group_2, conf_int, conf_level) {

  if (is.null(group_1) | is.null(group_2)) {
    group_1 <- levels(data$group)[1]
    group_2 <- levels(data$group)[2]
  }

  dat_1 <- data[data$group==group_1, ]
  dat_2 <- data[data$group==group_2, ]

  if (conf_int) {
    ratio_ci <- fieller_ratio_ci(a=dat_1$auc, b=dat_2$auc,
                                 a_se=dat_1$se, b_se=dat_2$se,
                                 conf_level=conf_level)

    out <- data.frame(to=dat_1$to,
                      ratio=ratio_ci$ratio,
                      ci_lower=ratio_ci$ci_lower,
                      ci_upper=ratio_ci$ci_upper)
    out$p_value <- fieller_p_val(a=dat_1$auc, b=dat_2$auc,
                                 a_se=dat_1$se, b_se=dat_2$se)
  } else {
    out <- data.frame(to=dat_1$to, ratio=dat_1$auc / dat_2$auc)
  }
  return(out)
}

## function to calculate the restricted mean survival time of each
## adjusted survival curve previously estimated using the adjustedsurv function
#' @export
adjusted_rmst <- function(adjsurv, to, from=0, conf_int=FALSE,
                          conf_level=0.95, interpolation="steps",
                          difference=FALSE, ratio=FALSE,
                          group_1=NULL, group_2=NULL) {

  check_inputs_adj_rmst(adjsurv=adjsurv, from=from, to=to, conf_int=conf_int,
                        difference=difference, ratio=ratio)

  # set to FALSE if it can't be done
  if (conf_int & is.null(adjsurv$boot_adjsurv)) {
    conf_int <- FALSE
  }

  out <- area_under_curve(adj=adjsurv, to=to, from=from,
                          conf_int=conf_int, conf_level=conf_level,
                          interpolation=interpolation)
  if (difference) {
    out <- auc_difference(data=out, group_1=group_1, group_2=group_2,
                          conf_int=conf_int, conf_level=conf_level)
  } else if (ratio) {
    out <- auc_ratio(data=out, group_1=group_1, group_2=group_2,
                     conf_int=conf_int, conf_level=conf_level)
  } else if (conf_int) {
    colnames(out) <- c("to", "group", "rmst", "se", "ci_lower", "ci_upper",
                       "n_boot")
  } else {
    colnames(out) <- c("to", "group", "rmst")
  }

  return(out)
}

## function to calculate the restricted mean time lost of each
## adjusted CIF previously estimated using the adjustedsurv/adjustedcif function
#' @export
adjusted_rmtl <- function(adj, to, from=0, conf_int=FALSE,
                          conf_level=0.95, interpolation="steps",
                          difference=FALSE, ratio=FALSE,
                          group_1=NULL, group_2=NULL) {

  check_inputs_adj_rmtl(adj=adj, from=from, to=to, conf_int=conf_int,
                        difference=difference, ratio=ratio)

  # set to FALSE if it can't be done
  if (conf_int & is.null(adj$boot_adjsurv) & is.null(adj$boot_adjcif)) {
    conf_int <- FALSE
  }

  # calculate area under curve
  out <- area_under_curve(adj=adj, to=to, from=from, conf_int=conf_int,
                          conf_level=conf_level,
                          interpolation=interpolation)

  # take area above curve instead
  if (inherits(adj, "adjustedsurv")) {

    full_area <- (to - from) * 1
    out$auc <- rep(full_area, each=length(unique(out$group))) - out$auc

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

  if (difference) {
    out <- auc_difference(data=out, group_1=group_1, group_2=group_2,
                          conf_int=conf_int, conf_level=conf_level)
  } else if (ratio) {
    out <- auc_ratio(data=out, group_1=group_1, group_2=group_2,
                     conf_int=conf_int, conf_level=conf_level)
  } else if (conf_int) {
    colnames(out) <- c("to", "group", "rmtl", "se", "ci_lower",
                       "ci_upper", "n_boot")
  } else {
    colnames(out) <- c("to", "group", "rmtl")
  }

  return(out)
}
