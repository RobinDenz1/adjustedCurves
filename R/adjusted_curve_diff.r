
## the function that actually does the estimation of both differences and
## ratios under the hood
adjusted_curve_diff.fit <- function(adj, group_1, group_2, times,
                                    conf_int, conf_level, use_boot,
                                    interpolation, type) {

  check_inputs_adj_diff(adj=adj, group_1=group_1, group_2=group_2,
                        conf_int=conf_int, use_boot=use_boot)

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

  # get difference / ratio + confidence intervals
  diff_curve <- difference_function(adj=plotdata, times=time_points,
                                    est=est_type, interpolation=interpolation,
                                    conf_int=conf_int, conf_level=conf_level,
                                    type=type, p_value=conf_int)

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

## calculate the curve of the difference including
## confidence intervals and p-values
#' @export
adjusted_curve_diff <- function(adj, group_1=NULL, group_2=NULL, times=NULL,
                                conf_int=FALSE, conf_level=0.95,
                                use_boot=FALSE, interpolation="steps") {

  out <- adjusted_curve_diff.fit(adj=adj, group_1=group_1, group_2=group_2,
                                 times=times, conf_int=conf_int,
                                 conf_level=conf_level, use_boot=use_boot,
                                 interpolation=interpolation, type="diff")
  return(out)
}

## calculate the ratio of two curves including
## confidence intervals and p-values
#' @export
adjusted_curve_ratio <- function(adj, group_1=NULL, group_2=NULL, times=NULL,
                                 conf_int=FALSE, conf_level=0.95,
                                 use_boot=FALSE, interpolation="steps") {

  out <- adjusted_curve_diff.fit(adj=adj, group_1=group_1, group_2=group_2,
                                 times=times, conf_int=conf_int,
                                 conf_level=conf_level, use_boot=use_boot,
                                 interpolation=interpolation, type="ratio")
  return(out)
}
