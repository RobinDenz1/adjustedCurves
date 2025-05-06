
## Adjustment based on a weighted average of stratified Kaplan-Meier estimates
## using the method by Cupples et al.
#' @export
surv_strat_cupples <- function(data, variable, ev_time, event,
                               conf_int=FALSE, conf_level=0.95, times,
                               adjust_vars, reference=NULL) {

  # devtools checks
  . <- .COVARS <- time <- treat_group <- surv <- count <- NULL

  # additional variables needed
  # NOTE: This code assumes that there is no column named .ALL or .COVARS
  #       and that there are no tabs in the column names
  data$.ALL <- interaction(data[, c(variable, adjust_vars)], sep="\t")
  if (is.null(reference)) {
    reference <- data
  }
  reference$.COVARS <- interaction(reference[, adjust_vars], sep="\t")

  # Kaplan-Meier survival curve for each possible strata at
  # every event time
  plotdata <- surv_km(data=data,
                      variable=".ALL",
                      ev_time=ev_time,
                      event=event,
                      times=times,
                      conf_int=FALSE)$plotdata

  # add indicator for treatment group and remove said treatment group from
  # the 'group' variable
  plotdata$treat_group <- sub("\t.*", "", plotdata$group)
  plotdata$group <- sub(".*?\t", "", plotdata$group)

  # add count of 'group' at baseline, ignoring treatment group
  count_dat <- reference %>%
    dplyr::group_by(.COVARS) %>%
    dplyr::summarise(count = dplyr::n())
  colnames(count_dat)[1] <- c("group")

  # merge together
  plotdata <- merge(count_dat, plotdata, by="group", all.y=TRUE)

  # take weighted in each treatment group at each t, over strata
  plotdata <- plotdata %>%
    dplyr::group_by(., time, treat_group) %>%
    dplyr::summarise(surv=stats::weighted.mean(x=surv, w=count, na.rm=FALSE),
                     .groups="drop_last")
  colnames(plotdata) <- c("time", "group", "surv")

  # keep same order in data.frame
  plotdata <- data.frame(time=plotdata$time,
                         surv=plotdata$surv,
                         group=plotdata$group)

  output <- list(plotdata=as.data.frame(plotdata))
  class(output) <- "adjustedsurv.method"

  return(output)
}

## Adjustment based on a weighted average of stratified Kaplan-Meier estimates
## using the method by Amato (1988)
#' @export
surv_strat_amato <- function(data, variable, ev_time, event,
                             conf_int=FALSE, conf_level=0.95,
                             times=NULL, adjust_vars, reference=NULL) {

  # silence checks
  . <- group <- time <- wdj <- wcj <- wrj <- delta_wdj <- delta_wd <- wr <- NULL
  times_input <- times

  # proportions in reference data
  if (is.null(reference)) {
    reference <- data
  }
  reference$.COVARS <- interaction(reference[, adjust_vars])
  Pjs <- prop.table(table(reference$.COVARS))

  # also calculate strata variable in data
  data$.COVARS <- interaction(data[, adjust_vars])

  levs <- levels(data[, variable])
  levs_adjust_var <- levels(data$.COVARS)
  out <- list()
  for (i in seq_len(length(levs))) {

    # data for treatment i
    dat_I <- data[data[, variable]==levs[i], ]
    times <- c(0, sort(unique(dat_I[, ev_time][dat_I[, event]==1])))

    for (j in seq_len(length(levs_adjust_var))) {

      # data for treatment i and only strata j
      dat_IJ <- dat_I[dat_I$.COVARS==levs_adjust_var[j], ]

      # weights for these individuals
      ajn <- nrow(dat_I) * Pjs[levs_adjust_var[j]] / nrow(dat_IJ)

      # people observed to fail
      Ndj <- vapply(times, FUN=function(x) {sum(dat_IJ[, ev_time] <= x &
                                                  dat_IJ[, event]==1)},
                    FUN.VALUE=numeric(1))
      delta_Ndj <- vapply(times, FUN=function(x) {sum(dat_IJ[, ev_time] == x &
                                                        dat_IJ[, event]==1)},
                          FUN.VALUE=numeric(1))

      # people censored
      Ncj <- vapply(times, FUN=function(x) {sum(dat_IJ[, ev_time] <= x &
                                                  dat_IJ[, event]==0)},
                    FUN.VALUE=numeric(1))
      # people at risk
      Nrj <- vapply(times, FUN=function(x) {sum(dat_IJ[, ev_time] >= x)},
                    FUN.VALUE=numeric(1))

      # put together
      temp <- data.frame(time=times, ajn=ajn[[1]],
                         Ndj=Ndj, Ncj=Ncj, Nrj=Nrj, delta_Ndj=delta_Ndj,
                         strata=levs_adjust_var[j], group=levs[i])
      out[[length(out)+1]] <- temp
    }
  }
  dat_stats <- dplyr::bind_rows(out)

  # calculate sums of weights
  dat_stats$wdj <- dat_stats$Ndj * dat_stats$ajn
  dat_stats$wcj <- dat_stats$Ncj * dat_stats$ajn
  dat_stats$wrj <- dat_stats$Nrj * dat_stats$ajn
  dat_stats$delta_wdj <- dat_stats$delta_Ndj * dat_stats$ajn

  # calculate survival probability
  plotdata <- dat_stats %>%
    dplyr::group_by(., group, time) %>%
    dplyr::summarise(wd=sum(wdj),
                     wc=sum(wcj),
                     wr=sum(wrj),
                     delta_wd=sum(delta_wdj),
                     .groups="drop_last") %>%
    dplyr::mutate(., surv=cumprod(1 - (delta_wd / wr)))

  # remove unnecessary variables
  plotdata$wd <- NULL
  plotdata$wc <- NULL
  plotdata$wr <- NULL
  plotdata$delta_wd <- NULL
  plotdata$group <- factor(plotdata$group, levels=levs)
  plotdata <- as.data.frame(plotdata)

  if (!is.null(times_input)) {
    plotdata <- specific_times(plotdata, times_input)
  }

  # keep same order in data.frame
  plotdata <- data.frame(time=plotdata$time,
                         surv=plotdata$surv,
                         group=factor(plotdata$group, levels=levs))

  output <- list(plotdata=plotdata,
                 Pjs=Pjs)
  class(output) <- "adjustedsurv.method"
  return(output)
}

## Adjustment based on a weighted average of stratified Kaplan-Meier estimates
## using the method by Gregory (1988) and Nieto & Coresh (1996)
# NOTE: Equations are due to Nieto & Coresh (1996) because while both
#       methods produce the same results when using the full data as reference,
#       only Nieto's formulation allows the calculation of confidence intervals.
#' @export
surv_strat_nieto <- function(data, variable, ev_time, event,
                             conf_int, conf_level=0.95,
                             times=NULL, adjust_vars) {

  # silence checks
  . <- time <- group <- frac <- est_var <- wji <- var_j <- NULL

  data$.COVARS <- interaction(data[, adjust_vars])
  times_input <- times

  # needed levels
  levs <- levels(data[, variable])
  levs_adjust_var <- levels(data$.COVARS)

  out <- list()
  for (i in seq_len(length(levs))) {

    dat_X <- data[data[, variable]==levs[i], ]

    # 1.)
    tj <- c(0, sort(unique(dat_X[, ev_time][dat_X[, event]==1])))

    for (j in seq_len(length(levs_adjust_var))) {

      # data at X, Z
      dat_XZ <- data[data[, variable]==levs[i] &
                       data$.COVARS==levs_adjust_var[j], ]

      # 2.)
      nxzj <- vapply(tj, function(x) {sum(dat_XZ[, ev_time]>=x)},
                     FUN.VALUE=numeric(1))
      axzj <- vapply(tj, function(x) {sum(dat_XZ[, ev_time]==x &
                                            dat_XZ[, event]==1)},
                     FUN.VALUE=numeric(1))
      # 3.)
      qxzj <- axzj / nxzj

      # 4.) but modified, calculating n at risk in strata overall instead
      #     of using the control group
      dat_Z <- data[data$.COVARS==levs_adjust_var[j], ]
      nz <- vapply(tj, function(x) {sum(dat_Z[, ev_time]>=x)},
                   FUN.VALUE=numeric(1))

      out[[length(out)+1]] <- data.frame(time=tj, nxzj=nxzj, axzj=axzj,
                                         qxzj=qxzj, nz=nz, group=levs[i],
                                         strata=levs_adjust_var[j])
    }
  }
  out <- dplyr::bind_rows(out)

  # appendix 1
  if (conf_int) {
    # calculate total nj (n at risk) in full data
    tj_overall <- c(0, sort(unique(data[, ev_time][data[, event]==1])))
    nj <- vapply(tj_overall, function(x) {sum(data[, ev_time]>=x)},
                 FUN.VALUE=numeric(1))
    dat_nj <- data.frame(time=tj_overall, nj=nj)

    # merge to previous output
    out <- merge(out, dat_nj, by="time", all.x=TRUE)

    # calculate wji
    out$wji <- out$nz / out$nj

    # calculate the sum needed at the end of equation 7
    out <- out %>%
      dplyr::group_by(., time, group) %>%
      dplyr::mutate(sum_wq=sum(wji * qxzj))

    # stratum specific variance
    out$var_j <- out$wji^2 * (((1 - out$qxzj)*out$qxzj) / out$nxzj) *
      (1 / (1 - out$sum_wq)^2)
  } else {
    out$var_j <- 0
  }

  # 5.) + 6.) and variance from appendix
  plotdata <- out %>%
    dplyr::group_by(., time, group) %>%
    dplyr::summarise(frac= (sum(qxzj * nz)) / sum(nz),
                     est_var=sum(var_j),
                     .groups="drop_last") %>%
    dplyr::group_by(., group) %>%
    dplyr::mutate(surv=cumprod(1 - frac),
                  se=sqrt(cumsum(est_var)))
  plotdata$frac <- NULL
  plotdata <- as.data.frame(plotdata)

  # equation 8 appendix
  if (conf_int) {
    surv_ci <- confint_surv(surv=log(plotdata$surv),
                            se=plotdata$se,
                            conf_level=conf_level,
                            conf_type="plain")
    plotdata$ci_lower <- exp(surv_ci$left)
    plotdata$ci_upper <- exp(surv_ci$right)
  } else {
    plotdata$se <- NULL
  }
  plotdata$est_var <- NULL

  if (!is.null(times_input)) {
    plotdata <- specific_times(plotdata, times_input)
  }

  plotdata_out <- data.frame(time=plotdata$time,
                             surv=plotdata$surv,
                             group=plotdata$group)
  plotdata_out$se <- plotdata$se
  plotdata_out$ci_lower <- plotdata$ci_lower
  plotdata_out$ci_upper <- plotdata$ci_upper

  output <- list(plotdata=plotdata)
  class(output) <- "adjustedsurv.method"

  return(output)
}
