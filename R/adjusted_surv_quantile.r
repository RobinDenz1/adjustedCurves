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

## get bootstrap standard error for differences or ratios when using
## multiply imputed data
#' @importFrom dplyr %>%
get_surv_quantile_boot_se.mi <- function(mids_analyses, p, interpolation) {

  # silcence devtools checks
  group <- se <- NULL

  # get standard error for all multiply imputed datasets
  out <- vector(mode="list", length=length(mids_analyses))
  for (i in seq_len(length(mids_analyses))) {
    boot_data_i <- mids_analyses[[i]]$boot_data
    out[[i]] <- get_surv_quantile_boot_se(boot_surv=boot_data_i,
                                          p=p, interpolation=interpolation)
  }
  out <- dplyr::bind_rows(out)

  # pool them using Rubins Rule
  boot_stats <- out %>%
    dplyr::group_by(p, group) %>%
    dplyr::summarise(se=mean(se, na.rm=TRUE),
                     .groups="drop_last")

  return(boot_stats)
}

## get bootstrap standard error for differences or ratios
#' @importFrom dplyr %>%
get_surv_quantile_boot_se <- function(boot_surv, p, interpolation) {

  # silence devtools checks
  q_surv <- group <- NULL

  # calculate survival time quantile for each bootstrap sample
  out <- vector(mode="list", length=max(boot_surv$boot))
  for (i in seq_len(max(boot_surv$boot))) {

    dat_i <- boot_surv[boot_surv$boot==i, ]

    q_boot_i <- get_surv_quantiles(plotdata=dat_i, p=p, conf_int=FALSE,
                                   interpolation=interpolation,
                                   conf_level=0.95, orig_conf_level=0.95,
                                   method="whatever", use_boot=FALSE)
    out[[i]] <- q_boot_i
  }
  out <- dplyr::bind_rows(out)

  # calculate standard error
  boot_stats <- out %>%
    dplyr::group_by(p, group) %>%
    dplyr::summarise(se=stats::sd(q_surv, na.rm=TRUE),
                     n_boot=sum(!is.na(q_surv)),
                     .groups="drop_last")

  return(boot_stats)
}

## core function to actually read off the survival time quantiles
#' @importFrom dplyr %>%
get_surv_quantiles <- function(plotdata, p, conf_int, conf_level,
                               interpolation, orig_conf_level,
                               method, use_boot, boot_data=NULL) {

  # silence devtools checks
  . <- time <- group <- surv <- NULL

  # remove NAs
  plotdata <- plotdata[!is.na(plotdata$surv),]
  levs <- levels(plotdata$group)

  # re-calculate confidence intervals to make sure that the
  # confidence levels match
  if (conf_int && (conf_level != orig_conf_level)) {

    if (method=="km") {
      stop("The 'conf_level' argument used in adjusted_surv_quantile()",
           " and the 'conf_level' argument used in adjustedsurv()",
           " must be set to the same value when using method='km'.",
           " Change 'conf_level' argument accordingly and run the code",
           " again for valid results.")
    }

    # using bootstrap percentiles
    if (use_boot) {
      boot_stats <- boot_data %>%
        dplyr::group_by(., time, group) %>%
        dplyr::summarise(ci_lower=stats::quantile(surv,
                                                  probs=(1-conf_level)/2,
                                                  na.rm=TRUE),
                         ci_upper=stats::quantile(surv,
                                                  probs=1-((1-conf_level)/2),
                                                  na.rm=TRUE),
                         .groups="drop_last")
      boot_stats <- as.data.frame(boot_stats)

      # order both data.frames just to be sure
      plotdata <- plotdata[order(plotdata$group, plotdata$time), ]
      boot_stats <- boot_stats[order(boot_stats$group,
                                     boot_stats$time), ]

      # put together
      boot_stats$surv <- plotdata$surv
      plotdata <- boot_stats

    # or using normal approximation
    } else {
      surv_ci <- confint_surv(surv=plotdata$surv,
                              se=plotdata$se,
                              conf_level=conf_level,
                              conf_type="plain")
      plotdata$ci_lower <- surv_ci$left
      plotdata$ci_upper <- surv_ci$right
    }
  }

  ## estimate the survival time quantiles by reading them off the survival
  ## curves directly
  out <- vector(mode="list", length=length(levs))
  for (i in seq_len(length(levs))) {
    temp_dat <- plotdata[plotdata$group==levs[i],]
    temp_dat$group <- NULL

    val <- read_from_fun(x=p, data=temp_dat, est="time", time="surv",
                         interpolation=interpolation)
    out_i <- data.frame(p=p, group=levs[i], q_surv=val)

    if (conf_int) {
      temp_dat <- temp_dat[!is.na(temp_dat$ci_lower) &
                             !is.na(temp_dat$ci_upper) ,]

      out_i$ci_lower <- read_from_fun(x=p, data=temp_dat,
                                      est="time", time="ci_lower",
                                      interpolation=interpolation)
      out_i$ci_upper <- read_from_fun(x=p, data=temp_dat,
                                      est="time", time="ci_upper",
                                      interpolation=interpolation)
    }
    out[[i]] <- out_i
  }
  out <- dplyr::bind_rows(out)

  return(out)
}

## get difference / ratio from survival time quantiles
get_surv_quantile_contrast <- function(plotdata, q_surv, contrast,
                                       group_1, group_2, p, conf_int,
                                       conf_level) {

  # define groups
  if (is.null(group_1) | is.null(group_2)) {
    levs <- levels(plotdata$group)
    group_1 <- levs[1]
    group_2 <- levs[2]
  }

  dat_1 <- q_surv[q_surv$group==group_1, ]
  dat_2 <- q_surv[q_surv$group==group_2, ]

  if (contrast=="diff") {
    out <- data.frame(p=p, diff=dat_1$q_surv - dat_2$q_surv)

    if (conf_int) {
      diff_se <- sqrt(dat_1$se^2 + dat_2$se^2)

      diff_ci <- confint_surv(surv=out$diff, se=diff_se,
                              conf_level=conf_level,
                              conf_type="plain")

      out$se <- diff_se
      out$ci_lower <- diff_ci$left
      out$ci_upper <- diff_ci$right

      t_stat <- (out$diff - 0) / out$se
      out$p_value <- 2 * stats::pnorm(abs(t_stat), lower.tail=FALSE)
    }

  } else if (contrast=="ratio") {
    out <- data.frame(p=p, ratio=dat_1$q_surv / dat_2$q_surv)

    if (conf_int) {
      ci_ratio <- fieller_ratio_ci(a=dat_1$q_surv, b=dat_2$q_surv,
                                   a_se=dat_1$se, b_se=dat_2$se,
                                   conf_level=conf_level)
      out$ci_lower <- ci_ratio$ci_lower
      out$ci_upper <- ci_ratio$ci_upper
      out$p_value <- fieller_p_val(a=dat_1$q_surv, b=dat_2$q_surv,
                                   a_se=dat_1$se, b_se=dat_2$se)
    }
  }

  return(out)
}

## calculate confounder adjusted survival time quantile
#' @export
adjusted_surv_quantile <- function(adjsurv, p=0.5, conf_int=FALSE,
                                   conf_level=adjsurv$conf_level,
                                   use_boot=FALSE, interpolation="steps",
                                   contrast="none", group_1=NULL,
                                   group_2=NULL) {

  check_inputs_surv_q(adjsurv=adjsurv, conf_int=conf_int, p=p,
                      use_boot=use_boot, interpolation=interpolation,
                      contrast=contrast, conf_level=conf_level,
                      group_1=group_1, group_2=group_2)

  if (use_boot) {
    plotdata <- adjsurv$boot_adj
  } else {
    plotdata <- adjsurv$adj
  }

  if (contrast %in% c("diff", "ratio")) {
    conf_int_main <- FALSE
  } else if (conf_int) {
    conf_int_main <- TRUE
  } else {
    conf_int_main <- FALSE
  }

  # estimate the quantiles
  out <- get_surv_quantiles(plotdata=plotdata, p=p, conf_int=conf_int_main,
                            interpolation=interpolation, conf_level=conf_level,
                            orig_conf_level=adjsurv$conf_level,
                            method=adjsurv$method, boot_data=adjsurv$boot_data,
                            use_boot=use_boot)

  if (contrast %in% c("diff", "ratio")) {
    if (conf_int) {

      if (is.null(adjsurv$mids_analyses)) {
        boot_se <- get_surv_quantile_boot_se(boot_surv=adjsurv$boot_data, p=p,
                                             interpolation=interpolation)
      } else {
        boot_se <- get_surv_quantile_boot_se.mi(
          mids_analyses=adjsurv$mids_analyses,
          interpolation=interpolation,
          p=p
        )
      }

      # put together
      out <- out[order(out$p, out$group), ]
      boot_se <- boot_se[order(boot_se$p, boot_se$group), ]
      out$se <- boot_se$se
    }

    out <- get_surv_quantile_contrast(plotdata=plotdata, q_surv=out, p=p,
                                      contrast=contrast,
                                      group_1=group_1, group_2=group_2,
                                      conf_int=conf_int,
                                      conf_level=conf_level)
  }
  return(out)
}
