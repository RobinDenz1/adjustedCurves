## Main function of the package. Is basically a wrapper around
## all other functions, offering additional high level stuff
# TODO: - might be nice to allow multicore execution for bootstrapping
#' @export
adjustedsurv <- function(data, variable, ev_time, event, method, sd=T,
                         times=NULL, alpha=0.05, bootstrap=F,
                         n_boot=500, na.rm=F, ...) {

  check_inputs_adjustedsurv(data=data, variable=variable,
                            ev_time=ev_time, event=event, method=method,
                            sd=sd, times=times, bootstrap=bootstrap,
                            n_boot=n_boot, na.rm=na.rm, ...)
  # get event specific times
  times_input <- times
  if (is.null(times)) {
    times <- sort(unique(data[, ev_time][data[, event]==1]))
  }

  levs <- unique(data[,variable])

  # get relevant surv_method function
  surv_fun <- get(paste0("surv_method_", method))

  # bootstrap the whole procedure, can be useful to get sd, p-values
  if (bootstrap) {

    boot_data_list <- list()
    boot_stats_list <- list()

    for (i in 1:n_boot) {
      indices <- sample(x=rownames(data), size=nrow(data), replace=T)
      boot_samp <- data[indices,]

      # if event specific times are used, use event specific times
      # in bootstrapping as well
      if (is.null(times_input)) {
        times_boot <- sort(unique(boot_samp[, ev_time][boot_samp[, event]==1]))
      } else {
        times_boot <- times
      }

      # update models/recalculate weights using bootstrap sample
      pass_args <- list(...)
      if (method %in% c("direct", "aiptw") &
          !inherits(pass_args$outcome_model, "formula")) {
        pass_args$outcome_model <- stats::update(pass_args$outcome_model,
                                                 data=boot_samp)
      }
      if (method %in% c("iptw_km", "iptw_cox", "iptw_pseudo", "aiptw",
                        "aiptw_pseudo")) {
        pass_args$treatment_model <- stats::update(pass_args$treatment_model,
                                                   data=boot_samp, trace=F)
      }

      # call surv_method with correct arguments
      args <- list(data=boot_samp, variable=variable, ev_time=ev_time,
                   event=event, sd=F, times=times, na.rm=na.rm)
      args <- c(args, pass_args)

      adjsurv_boot <- R.utils::doCall(surv_fun, args=args)
      adjsurv_boot$boot <- i
      boot_data_list[[i]] <- adjsurv_boot

      # read from resulting step function at all t in times
      boot_surv <- vector(mode="list", length=length(levs))
      for (j in 1:length(levs)) {

        if (method %in% c("iptw_km", "iptw_cox")) {
          times <- unique(data[,ev_time][data[,variable]==levs[j]])
        }

        surv_boot <- sapply(times, read_from_step_function,
                            step_data=adjsurv_boot[adjsurv_boot$group==levs[j],])

        dat_temp <- data.frame(time=times,
                               surv_b=surv_boot,
                               group=levs[j],
                               boot=i)
        boot_stats_list[[length(boot_stats_list)+1]] <- dat_temp
      }
    }
    boot_data <- as.data.frame(dplyr::bind_rows(boot_data_list))
    boot_data_same_t <- as.data.frame(dplyr::bind_rows(boot_stats_list))

    # calculate some statistics
    boot_stats <- boot_data_same_t %>%
      dplyr::group_by(., time, group) %>%
      dplyr::summarise(surv=mean(surv_b, na.rm=na.rm),
                       sd=stats::sd(surv_b, na.rm=na.rm),
                       ci_lower=stats::quantile(surv_b, probs=alpha/2,
                                                na.rm=na.rm),
                       ci_upper=stats::quantile(surv_b, probs=1-(alpha/2),
                                                na.rm=na.rm),
                       .groups="drop_last")
  }

  # core of the function
  args <- list(data=data, variable=variable, ev_time=ev_time,
               event=event, sd=sd, times=times, na.rm=na.rm, ...)
  plotdata <- R.utils::doCall(surv_fun, args=args)

  # calculate point-wise confidence intervals
  if (sd) {
    plotdata$ci_lower <- confint_surv(surv=plotdata$surv, sd=plotdata$sd,
                                      n=nrow(data), alpha=alpha,
                                      conf_type="plain")$left
    plotdata$ci_upper <- confint_surv(surv=plotdata$surv, sd=plotdata$sd,
                                      n=nrow(data), alpha=alpha,
                                      conf_type="plain")$right
  }

  out <- list(adjsurv=plotdata,
              method=method,
              categorical=ifelse(length(unique(data[,variable]))>2, T, F),
              call=match.call())

  if (bootstrap) {
    out$boot_data <- boot_data
    out$boot_data_same_t <- boot_data_same_t
    out$boot_adjsurv <- boot_stats
  }

  class(out) <- "adjustedsurv"
  return(out)
}

## plot the survival curves
# TODO: maybe need to recalculate confidence intervals when using iso_reg?
#' @importFrom rlang .data
#' @export
plot.adjustedsurv <- function(adjsurv, draw_ci=T, max_t=Inf,
                              iso_reg=F, force_bounds=F, use_boot=F,
                              color=T, linetype=F, facet=F,
                              line_size=1, xlab="Time",
                              ylab="Adjusted Survival Probability",
                              title=NULL, legend_title="Group",
                              legend_position="right",
                              ylim=NULL, custom_colors=NULL,
                              custom_linetypes=NULL,
                              ci_draw_alpha=0.4) {

  if (use_boot & is.null(adjsurv$boot_adjsurv)) {
    stop("Cannot use bootstrapped estimates as they were not estimated.",
         " Need bootstrap=TRUE in adjustedsurv() call.")
  } else if (use_boot) {
    plotdata <- adjsurv$boot_adjsurv
  } else {
    plotdata <- adjsurv$adjsurv
  }
  plotdata$group <- factor(plotdata$group)

  # shortcut to only show curves up to a certain time
  plotdata <- plotdata[which(plotdata$time <= max_t),]

  # in some methods estimates can be outside the 0, 1 bounds,
  # if specified set those to 0 or 1 respectively
  if (force_bounds) {
    plotdata <- within(plotdata, {
      surv <- ifelse(surv < 0, 0, surv);
      surv <- ifelse(surv > 1, 1, surv);
    })
  }

  # apply isotonic regression if specified
  if (iso_reg) {
    for (lev in levels(plotdata$group)) {
      surv <- plotdata$surv[plotdata$group==lev]
      new <- rev(stats::isoreg(rev(surv))$yf)
      plotdata$surv[plotdata$group==lev] <- new
    }
  }

  mapping <- ggplot2::aes(x=.data$time, y=.data$surv, color=.data$group,
                          linetype=.data$group)

  if (!linetype) {
    mapping$linetype <- NULL
  }
  if (!color) {
    mapping$colour <- NULL
  }

  p <- ggplot2::ggplot(plotdata, mapping) +
    ggplot2::geom_step(size=line_size) +
    ggplot2::theme_bw() +
    ggplot2::labs(x=xlab, y=ylab, color=legend_title,
                  linetype=legend_title, fill=legend_title) +
    ggplot2::theme(legend.position=legend_position)

  if (facet) {
    p <- p + ggplot2::facet_wrap(~group)
  }
  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }
  if (!is.null(ylim)) {
    p <- p + ggplot2::ylim(ylim)
  }
  if (!is.null(custom_colors)) {
    p <- p + ggplot2::scale_colour_manual(values=custom_colors)
  }
  if (!is.null(custom_linetypes)) {
    p <- p + ggplot2::scale_linetype_manual(values=custom_linetypes)
  }
  if (!is.null(custom_colors)) {
    p <- p + ggplot2::scale_fill_manual(values=custom_colors)
  }

  if (draw_ci & !"sd" %in% colnames(plotdata)) {
    warning("Cannot draw confidence intervals. Need 'sd=TRUE' in",
            " 'adjustedsurv()' call. Alternatively, use the bootstrap estimates.")
  }

  if (draw_ci & "sd" %in% colnames(plotdata)) {
    p <- p + pammtools::geom_stepribbon(ggplot2::aes(ymin=.data$ci_lower,
                                            ymax=.data$ci_upper,
                                            fill=.data$group,
                                            x=.data$time,
                                            y=.data$surv),
                                        alpha=ci_draw_alpha, inherit.aes=F)
  }
  return(p)
}

## Function to calculate some statistics for the survival curves
# TODO:
# - only allow estimation if both curves were estimated up to "to"
# - should work with pairwise comparisons
# - allow multicore / parallel processing
# - maybe write a plot method for this function
#' @export
adjustedsurv_test <- function(adjsurv, from=0, to=Inf) {

  if (!inherits(adjsurv, "adjustedsurv")) {
    stop("'adjsurv' must be an 'adjustedsurv' object, created using the ",
         "adjustedsurv function.")
  } else if (is.null(adjsurv$boot_data)) {
    stop("Can only perform a significance test if bootstrapping was ",
         "performed (bootstrap=TRUE in adjsutedsurv() call).")
  } else if (adjsurv$categorical) {
    stop("This function currently only supports a test of two survival curves.")
  }

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
                                          max_t=to)
    curve_list[[i]] <- surv_diff

    # integral of that curve
    diff_integral <- exact_stepfun_integral(surv_diff)

    stats_vec[i]<- diff_integral
  }

  diff_curves <- as.data.frame(dplyr::bind_rows(curve_list))

  ## use those bootstrapped integrals for the calculation of
  ## a p-value, by shifting the observed distribution to 0 and
  ## comparing it to the actually observed value

  # actually observed values
  times <- sort(unique(adjsurv$adjsurv$time))
  observed_diff_curve <- exact_stepfun_difference(adjsurv=adjsurv$adjsurv,
                                                  times=times,
                                                  max_t=to)
  observed_diff_integral <- exact_stepfun_integral(observed_diff_curve)

  # shit bootstrap distribution
  diff_under_H0 <- stats_vec - mean(stats_vec)
  p_value <- mean(abs(diff_under_H0) > abs(observed_diff_integral))

  # put together output
  out <- list(diff_curves=diff_curves,
              diff_integrals=stats_vec,
              observed_diff_curve=observed_diff_curve,
              observed_diff_integral=observed_diff_integral,
              p_value=p_value,
              n_boot=length(stats_vec),
              from=from,
              to=to)
  class(out) <- "adjustedsurv_test"

  return(out)
}

## print method for adjustedsurv_test
#' @export
print.adjustedsurv_test <- function(adjtest) {

  cat("##################################################################\n")
  cat("Pepe-Flemming Test of Equality of Two Adjusted Survival Curves \n")
  cat("------------------------------------------------------------------")
  cat("\n")
  cat("The equality was tested for the time interval: ", adjtest$from, " to ",
      adjtest$to, "\n")
  cat("Observed Integral of the difference: ", adjtest$observed_diff_integral,
      "\n")
  cat("Bootstrap standard deviation: ", stats::sd(adjtest$diff_integrals), "\n")
  cat("P-Value: ", adjtest$p_value, "\n\n")
  cat("Calculated using ", adjtest$n_boot, " bootstrap replications.\n")
  cat("##################################################################\n")
}

## function to calculate the restricted mean survival time of each
## adjusted survival curve previously estimated using the adjustedsurv function
# TODO: - allow multicore bootstrapping
#       - allow confidence interval calculation for rmsts
#' @export
adjusted_rmst <- function(adjsurv, from=0, to=Inf, use_boot=F) {

  check_inputs_adj_rmst(adjsurv=adjsurv, from=from, to=to, use_boot=use_boot)

  if (use_boot) {

    n_boot <- max(adjsurv$boot_data$boot)
    booted_areas <- vector(mode="list", length=n_boot)
    booted_rmsts <- vector(mode="list", length=n_boot)

    for (i in 1:n_boot) {

      # select one bootstrap data set each
      boot_dat <- adjsurv$boot_data[adjsurv$boot_data$boot==i,]

      # create fake adjustedsurv object
      fake_adjsurv <- list(adjsurv=boot_dat)
      class(fake_adjsurv) <- "adjustedsurv"

      # recursion call
      adj_rmst <- adjusted_rmst(fake_adjsurv, from=from, to=to, use_boot=F)

      booted_areas[[i]] <- adj_rmst$areas
      booted_rmsts[[i]] <- adj_rmst$rmsts
    }
    booted_areas <- dplyr::bind_rows(booted_areas)
    booted_rmsts <- dplyr::bind_rows(booted_rmsts)
  }

  levs <- unique(adjsurv$adjsurv$group)

  rmsts <- vector(mode="numeric", length=length(levs))
  areas <- vector(mode="numeric", length=length(levs))
  for (i in 1:length(levs)) {

    surv_dat <- adjsurv$adjsurv[adjsurv$adjsurv$group==levs[i],]
    surv_dat$group <- NULL
    surv_dat$sd <- NULL
    surv_dat$ci_lower <- NULL
    surv_dat$ci_upper <- NULL
    surv_dat$boot <- NULL

    # constrain step function end
    if (is.finite(to)) {
      latest <- read_from_step_function(to, step_data=surv_dat)
      surv_dat <- surv_dat[surv_dat$time <= to,]

      if (!to %in% surv_dat$time) {
        surv_dat <- rbind(surv_dat, data.frame(time=to, surv=latest))
      }

    }

    # constrain step function beginning
    if (from != 0) {
      earliest <- read_from_step_function(from, step_data=surv_dat)
      surv_dat <- surv_dat[surv_dat$time >= from,]

      if (!from %in% surv_dat$time) {
        surv_dat <- rbind(data.frame(time=from, surv=earliest), surv_dat)
      }

    } else {
      if (!0 %in% surv_dat$time) {
        surv-dat <- rbind(data.frame(time=0, surv=1), surv_dat)
      }
    }

    rmst <- exact_stepfun_integral(surv_dat)
    areas[i] <- rmst

    rmst <- 1/max(surv_dat$time) * rmst
    rmsts[i] <- rmst
  }
  names(rmsts) <- levs
  names(areas) <- levs

  out <- list(areas=areas,
              rmsts=rmsts,
              from=from,
              to=to)
  class(out) <- "adjusted_rmst"

  if (use_boot) {
    out$n_boot <- max(adjsurv$boot_data$boot)
    out$booted_areas <- booted_areas
    out$booted_rmsts <- booted_rmsts
    out$areas_sd <- apply(booted_areas, 2, sd)
    out$rmsts_sd <- apply(booted_rmsts, 2, sd)
  }

  return(out)
}

## print method for adjusted_rmst function
#' @export
print.adjusted_rmst <- function(adj_rmst, digits=5) {

  rmsts_str <- paste(names(adj_rmst$rmsts), round(adj_rmst$rmsts, digits),
                     sep="=", collapse="  ")
  areas_str <- paste(names(adj_rmst$areas), round(adj_rmst$areas, digits),
                     sep="=", collapse="  ")

  cat("Confounder-Adjusted Restricted Mean Survival Time\n")
  cat("\n")
  cat("Using the interval:", adj_rmst$from, "to", adj_rmst$to, "\n")
  cat("RMSTS: ", rmsts_str, "\n")
  cat("Areas under the curves (AUC): ", areas_str, "\n")

  if (!is.null(adj_rmst$booted_areas)) {
    rmsts_sd_str <- paste(names(adj_rmst$rmsts_sd), round(adj_rmst$rmsts_sd, digits),
                          sep="=", collapse="  ")
    areas_sd_str <- paste(names(adj_rmst$areas_sd), round(adj_rmst$areas_sd, digits),
                          sep="=", collapse="  ")

    cat("RMSTS Standard Deviation: ", rmsts_sd_str, "\n")
    cat("AUC Standard Deviation: ", areas_sd_str, "\n")
    cat("\n")
    cat("SD estimated using", adj_rmst$n_boot, "bootstrap replications.")
  }

}

## function to simulate confounded survival data
#' @export
sim_confounded_surv <- function(n=500, lcovars=NULL, outcome_betas=NULL,
                                group_beta=-1, surv_dist="weibull",
                                gamma=1.8, lambda=2, treatment_betas=NULL,
                                intercept=-0.5, gtol=0.001,
                                cens_fun=function(n){stats::rweibull(n, 1, 2)},
                                cens_args=list(), max_t=Inf) {

  check_inputs_sim_fun(n=n, lcovars=lcovars, outcome_betas=outcome_betas,
                       surv_dist=surv_dist, gamma=gamma, lambda=lambda,
                       treatment_betas=treatment_betas, group_beta=group_beta,
                       intercept=intercept, gtol=gtol, cens_fun=cens_fun,
                       cens_args=cens_args, max_t=max_t)

  # set defaults to parameters used in authors simulation study
  if (is.null(lcovars)) {
    lcovars <- list(x1=c("rbinom", 1, 0.5),
                    x2=c("rbinom", 1, 0.5),
                    x3=c("rbinom", 1, 0.5),
                    x4=c("rnorm", 0, 1, -2, 2),
                    x5=c("rnorm", 0, 1, -2, 2),
                    x6=c("rnorm", 0, 1, -2, 2))
  }
  if (is.null(outcome_betas)) {
    outcome_betas <- c(x1=log(1.8), x2=log(1.3),
                       x3=0, x4=log(1.8),
                       x5=log(1.3), x6=0)
  }
  if (is.null(treatment_betas)) {
    treatment_betas <- c(x1=0, x2=log(2), x3=log(1.5),
                         x4=0, x5=log(1.5), x6=log(2))
  }

  # draw i.i.d. random variables specified in lcovars
  covars <- data.frame(id=1:n)
  for (col in names(lcovars)) {
    if (lcovars[[col]][1]=="rnorm") {
      covars[,col] <- stats::rnorm(n, as.numeric(lcovars[[col]][2]),
                            as.numeric(lcovars[[col]][3]))
    } else if (lcovars[[col]][1]=="runif") {
      covars[,col] <- stats::runif(n, as.numeric(lcovars[[col]][2]),
                            as.numeric(lcovars[[col]][3]))
    } else if (lcovars[[col]][1]=="rbinom") {
      covars[,col] <- stats::rbinom(n, as.numeric(lcovars[[col]][2]),
                             as.numeric(lcovars[[col]][3]))
    }
  }
  covars$id <- NULL

  # assign binary treatment using logistic regression
  group_p <- intercept + rowSums(treatment_betas * covars[,names(treatment_betas)])
  group_p <- 1/(1 + exp(-group_p))

  # in order to keep the positivity assumption,
  # values of 0 and 1 can not be tolerated
  group_p[group_p < gtol] <- gtol
  group_p[group_p > (1 - gtol)] <- 1 - gtol

  covars$group <- stats::rbinom(n=n, size=1, prob=group_p)

  # generate survival times
  covars$time <- apply(X=covars, MARGIN=1, FUN=sim_surv_time,
                       betas=c(outcome_betas, group=group_beta),
                       dist=surv_dist, lambda=lambda, gamma=gamma)
  covars$event <- 1

  # introduce random censoring if specified
  if (!is.null(cens_fun)) {
    cens_time <- do.call(cens_fun, args=c(n=n, cens_args))

    covars <- within(covars, {
      event <- ifelse(time < cens_time, 1, 0);
      time <- ifelse(time < cens_time, time, cens_time);
    })
  }
  # also add administrative censoring if specified
  covars <- within(covars, {
    event <- ifelse(time <= max_t, event, 0);
    time <- ifelse(time <= max_t, time, max_t);
  })

  return(covars)
}
