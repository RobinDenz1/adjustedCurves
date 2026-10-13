
## the function that actually does the estimation of both differences and
## ratios under the hood
adjusted_curve_diff.fit <- function(adj, group_1, group_2, times,
                                    conf_int, conf_level, use_boot,
                                    interpolation, type, paired_boot) {

  check_inputs_adj_diff(adj=adj, group_1=group_1, group_2=group_2,
                        conf_int=conf_int, use_boot=use_boot,
                        paired_boot=paired_boot)

  # silence devtools::check()
  . <- time <- est <- NULL

  # get plotdata out of adj
  if (use_boot) {
    plotdata <- adj$boot_adj
  } else {
    plotdata <- adj$adj
  }

  if (inherits(adj, "adjustedsurv")) {
    est_type <- "surv"
  } else if (inherits(adj, "adjustedcif")) {
    est_type <- "cif"
  }

  # set groups
  if (is.null(group_1) | is.null(group_2)) {
    group_1 <- levels(plotdata$group)[1]
    group_2 <- levels(plotdata$group)[2]
  }

  # keep only data with the groups of interest
  plotdata <- plotdata[plotdata$group==group_1 | plotdata$group==group_2,]
  plotdata$group <- factor(plotdata$group, levels=c(group_1, group_2))

  # set default time points
  if (is.null(times)) {
    time_points <- sort(unique(plotdata$time))
  } else {
    time_points <- times
  }

  # use the (co-)variances of the paired bootstrap estimates, so that the
  # two curves do not have to be assumed independent
  boot_cov <- NULL
  if (conf_int && use_boot && paired_boot) {
    boot_cov <- paired_boot_cov(adj=adj, group_1=group_1, group_2=group_2,
                                times=time_points, est=est_type,
                                interpolation=interpolation)
  }

  # get difference / ratio + confidence intervals
  diff_curve <- difference_function(adj=plotdata, times=time_points,
                                    est=est_type, interpolation=interpolation,
                                    conf_int=conf_int, conf_level=conf_level,
                                    type=type, p_value=conf_int,
                                    boot_cov=boot_cov)

  # rename columns in output
  if (conf_int & type=="diff") {
    colnames(diff_curve) <- c("time", type, "se", "ci_lower", "ci_upper",
                              "p_value")
  } else if (conf_int & type=="ratio") {
    colnames(diff_curve) <- c("time", type, "ci_lower", "ci_upper",
                              "p_value")
  } else {
    colnames(diff_curve) <- c("time", type)
  }

  return(diff_curve)
}

## estimate the variances of the two curves and their covariance at some
## points in time from the bootstrap replications, in which both curves
## were estimated using the same bootstrap sample
paired_boot_cov <- function(adj, group_1, group_2, times, est,
                            interpolation) {

  # read one curve at all points in time, NA if it can't be read
  read_curve <- function(data) {
    rep_len(read_from_fun(x=times, data=data, est=est,
                          interpolation=interpolation), length(times))
  }

  # variance of x, variance of y and covariance, using complete pairs only
  cov_pairs <- function(x, y) {
    ok <- !is.na(x) & !is.na(y)
    if (sum(ok) < 2) {
      return(rep(NA_real_, 3))
    }
    c(stats::var(x[ok]), stats::var(y[ok]), stats::cov(x[ok], y[ok]))
  }

  # multiple imputation: pool the (co-)variances using Rubin's rule
  if (!is.null(adj$mids_analyses)) {

    m <- length(adj$mids_analyses)
    imp_covs <- lapply(adj$mids_analyses, paired_boot_cov, group_1=group_1,
                       group_2=group_2, times=times, est=est,
                       interpolation=interpolation)
    imp_est <- function(group) {
      matrix(vapply(adj$mids_analyses, function(x) {
        read_curve(x$boot_adj[x$boot_adj$group==group, ])
      }, numeric(length(times))), nrow=length(times))
    }
    est_1 <- imp_est(group_1)
    est_2 <- imp_est(group_2)

    within <- vapply(c("var_0", "var_1", "cov_01"), function(col) {
      rowMeans(matrix(vapply(imp_covs, function(x) x[, col],
                             numeric(length(times))), nrow=length(times)),
               na.rm=TRUE)
    }, numeric(length(times)))
    between <- t(vapply(seq_along(times), function(i) {
      cov_pairs(est_1[i, ], est_2[i, ])
    }, numeric(3)))
    total <- matrix(within, nrow=length(times)) + between * (1 + 1/m)

  } else {

    boot_data <- adj$boot_data
    boot_ids <- unique(boot_data$boot)
    boot_curves <- function(group) {
      dat <- boot_data[boot_data$group==group, ]
      dats <- split(dat, factor(dat$boot, levels=boot_ids))
      matrix(vapply(dats, read_curve, numeric(length(times))),
             nrow=length(times))
    }
    curves_1 <- boot_curves(group_1)
    curves_2 <- boot_curves(group_2)

    total <- t(vapply(seq_along(times), function(i) {
      cov_pairs(curves_1[i, ], curves_2[i, ])
    }, numeric(3)))
  }

  out <- data.frame(time=times, var_0=total[, 1], var_1=total[, 2],
                    cov_01=total[, 3])
  return(out)
}

## calculate the curve of the difference including
## confidence intervals and p-values
#' @export
adjusted_curve_diff <- function(adj, group_1=NULL, group_2=NULL, times=NULL,
                                conf_int=FALSE, conf_level=0.95,
                                use_boot=FALSE, interpolation="steps",
                                paired_boot=TRUE) {

  out <- adjusted_curve_diff.fit(adj=adj, group_1=group_1, group_2=group_2,
                                 times=times, conf_int=conf_int,
                                 conf_level=conf_level, use_boot=use_boot,
                                 interpolation=interpolation, type="diff",
                                 paired_boot=paired_boot)
  return(out)
}

## calculate the ratio of two curves including
## confidence intervals and p-values
#' @export
adjusted_curve_ratio <- function(adj, group_1=NULL, group_2=NULL, times=NULL,
                                 conf_int=FALSE, conf_level=0.95,
                                 use_boot=FALSE, interpolation="steps",
                                 paired_boot=TRUE) {

  out <- adjusted_curve_diff.fit(adj=adj, group_1=group_1, group_2=group_2,
                                 times=times, conf_int=conf_int,
                                 conf_level=conf_level, use_boot=use_boot,
                                 interpolation=interpolation, type="ratio",
                                 paired_boot=paired_boot)
  return(out)
}
