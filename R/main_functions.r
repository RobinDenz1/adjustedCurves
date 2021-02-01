###################### Main, higher level functions ############################

#setwd("C:/Users/Robin Denz/Documents/Arbeit/R-Paket_adjustedcurves")

#source("helper_functions.r")

# BIG TODO:
# - make times object work as intended
# - make it all work with factors
#     - matching, pseudo_dr, tmle, ostmle
# - allow categorical if appropriate
# - calculate simultaneous confidence bands (not sure how yet)
# - allow bootstrap CI for all methods
# - check if the CIs i have now make sense with simulation

# Needed functions:
# - function to calculate area between curves
# - maybe even bootstrap based significance testing
# - LATER: CIF methods with similar functionality
# - need unit tests (testthat maybe?)

## example values for testing
# censoring function
#custom_cens <- function(n) {
#  cens_times <- rweibull(n, 1, 2)
#  cens_times <- round(cens_times*1000) + 1
#  return(cens_times)
#}
#
#data <- sample_from_pop(pop, n=500, cens_fun=custom_cens, max_t=1826)
#data$group <- factor(data$group)
#variable <- "group"
#ev_time <- "time"
#event <- "cause"
#sd = T
#treat_model <- glm(group ~ x2 + x4 + x5, data=data, family="binomial")
#ps_score <- treat_model$fitted.values
#outcome_model <- coxph(Surv(time, cause) ~ ., data=data, x=T)
#outcome_vars <- c("x1", "x2", "x3", "x4", "x5", "x6")
#SL.trt <- c("SL.glm")
#SL.ftime <- c("SL.glm")
#SL.ctime <- c("SL.glm")

#adjsurv <- adjusted_surv(data=data, variable=variable, ev_time=ev_time,
#                         event=event, method="direct", sd=T,
#                         outcome_model=outcome_model)

## Main function
#' @export
adjusted_surv <- function(data, variable, ev_time, event, method, sd=T,
                          times="events_overall", alpha=0.05, ...) {

  check_inputs_adjusted_surv(data=data, variable=variable,
                             ev_time=ev_time, event=event, method=method,
                             sd=sd, times=times, ...)

  times <- get_times(data=data, variable=variable, ev_time=ev_time,
                     event=event, type=times, custom=NULL)

  surv_fun <- get(paste0("surv_method_", method))

  if (method %in% c("km", "matching")) {
    args <- list(data=data, variable=variable, ev_time=ev_time,
                 event=event, sd=sd, ...)
  } else {
    args <- list(data=data, variable=variable, ev_time=ev_time,
                 event=event, sd=sd, times=times, ...)
  }
  plotdata <- do.call(surv_fun, args=args)

  # calculate point-wise confidence intervals
  if (sd) {
    plotdata <- within(plotdata, {
      ci_lower <- confint_surv(surv=surv, sd=sd, n=nrow(data),
                               alpha=alpha)$left;
      ci_upper <- confint_surv(surv=surv, sd=sd, n=nrow(data),
                               alpha=alpha)$right;
    })
  }

  out <- list(adjsurv=plotdata,
              method=method,
              call=match.call())
  class(out) <- "adjustedsurv"
  return(out)
}

## plot the survival curves
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_step theme_bw theme facet_wrap
#' ggtitle ylim scale_colour_manual scale_linetype_manual
#' @export
plot.adjustedsurv <- function(adjsurv, draw_ci=T, max_t=Inf,
                              isotonic_regression=F, force_bounds=F,
                              color=T, linetype=F, facet=F,
                              line_size=1, xlab="Time",
                              ylab="Adjusted Survival Probability",
                              title=NULL, legend_title="Group",
                              legend_position="right",
                              ylim=NULL, custom_colors=NULL,
                              custom_linetypes=NULL,
                              ci_draw_alpha=0.4) {

  plotdata <- adjsurv$adjsurv
  plotdata$group <- factor(plotdata$group)

  # shortcut to only show curves up to a certain time
  plotdata <- subset(plotdata, time <= max_t)

  # in some methods estimates can be outside the 0, 1 bounds,
  # if specified set those to 0 or 1 respectively
  if (force_bounds) {
    plotdata <- within(plotdata, {
      surv <- ifelse(surv < 0, 0, surv);
      surv <- ifelse(surv > 1, 1, surv);
    })
  }

  # apply isotonic regression if specified
  # TODO: maybe need to recalculate confidence intervals
  if (isotonic_regression) {
    for (lev in levels(plotdata$group)) {
      surv <- plotdata$surv[plotdata$group==lev]
      new <- rev(stats::isoreg(rev(surv))$yf)
      plotdata$surv[plotdata$group==lev] <- new
    }
  }

  mapping <- aes(x=.data$time, y=.data$surv, color=.data$group,
                 linetype=.data$group)

  if (!linetype) {
    mapping$linetype <- NULL
  }
  if (!color) {
    mapping$colour <- NULL
  }

  p <- ggplot(plotdata, mapping) +
    geom_step(size=line_size) +
    theme_bw() +
    labs(x=xlab, y=ylab, color=legend_title,
         linetype=legend_title, fill=legend_title) +
    theme(legend.position=legend_position)

  if (facet) {
    p <- p + facet_wrap(~group)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  if (!is.null(ylim)) {
    p <- p + ylim(ylim)
  }
  if (!is.null(custom_colors)) {
    p <- p + scale_colour_manual(values=custom_colors)
  }
  if (!is.null(custom_linetypes)) {
    p <- p + scale_linetype_manual(values=custom_linetypes)
  }
  if (!is.null(custom_colors)) {
    p <- p + scale_fill_manual(values=custom_colors)
  }

  if (draw_ci & !"sd" %in% colnames(plotdata)) {
    warning("Cannot draw confidence intervals. Need 'sd=TRUE' in",
            " 'adjusted_surv()' call.")
  }

  if (draw_ci & "sd" %in% colnames(plotdata)) {
    p <- p + pammtools::geom_stepribbon(aes(ymin=.data$ci_lower,
                                            ymax=.data$ci_upper,
                                            fill=.data$group,
                                            x=.data$time,
                                            y=.data$surv),
                                        alpha=ci_draw_alpha, inherit.aes=F)
  }
  return(p)
}
