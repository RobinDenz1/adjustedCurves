## Main function of the package. Is basically a wrapper around
## all other functions, offering additional high level stuff
#' @importFrom dplyr %>%
#' @importFrom doRNG %dorng%
#' @importFrom foreach %dopar%
#' @export
adjustedsurv <- function(data, variable, ev_time, event, method, conf_int=F,
                         conf_level=0.95, times=NULL, bootstrap=F,
                         n_boot=500, n_cores=1, ...) {

  check_inputs_adjustedsurv(data=data, variable=variable,
                            ev_time=ev_time, event=event, method=method,
                            conf_int=conf_int, conf_level=conf_level,
                            times=times, bootstrap=bootstrap,
                            n_boot=n_boot, ...)

  # define those to remove Notes in devtools::check()
  . <- i <- time <- group <- surv_b <- NULL

  # get event specific times
  times_input <- times
  if (is.null(times)) {
    times <- sort(unique(data[, ev_time][data[, event]==1]))

    # add zero if not already in there
    if (!0 %in% times) {
      times <- c(0, times)
    }
  }

  levs <- unique(data[,variable])

  # get relevant surv_method function
  surv_fun <- get(paste0("surv_", method))

  # bootstrap the whole procedure, can be useful to get sd, p-values
  if (bootstrap) {

    if (n_cores > 1) {
      # needed packages for parallel processing
      requireNamespace("parallel")
      requireNamespace("doRNG")
      requireNamespace("doParallel")
      requireNamespace("foreach")

      # initialize clusters
      cl <- parallel::makeCluster(n_cores, outfile="")
      doParallel::registerDoParallel(cl)
      pkgs <- c("adjustedCurves", "survival")
      export_objs <- c("get_iptw_weights", "read_from_step_function",
                       "multi_result_class", "adjustedsurv_boot")

      boot_out <- foreach::foreach(i=1:n_boot, .packages=pkgs,
                                   .export=export_objs) %dorng% {

      adjustedsurv_boot(data=data, variable=variable, ev_time=ev_time,
                        event=event, method=method, times_input=times_input,
                        times=times, i=i, surv_fun=surv_fun,
                        levs=levs, ...)
      }
      parallel::stopCluster(cl)

    } else {

      boot_out <- vector(mode="list", length=n_boot)
      for (i in 1:n_boot) {
        boot_out[[i]] <- adjustedsurv_boot(data=data, variable=variable,
                                           ev_time=ev_time, event=event,
                                           method=method, times_input=times_input,
                                           times=times, i=i,
                                           surv_fun=surv_fun, levs=levs, ...)
      }
    }

    # transform into data.frames
    boot_data <- lapply(boot_out, function(x) x$boot_data)
    boot_data_same_t <- lapply(boot_out, function(x) x$boot_data_same_t)

    boot_data <- as.data.frame(dplyr::bind_rows(boot_data))
    boot_data_same_t <- as.data.frame(dplyr::bind_rows(boot_data_same_t))

    # calculate some statistics
    boot_stats <- boot_data_same_t %>%
      dplyr::group_by(., time, group) %>%
      dplyr::summarise(surv=mean(surv_b, na.rm=T),
                       sd=stats::sd(surv_b, na.rm=T),
                       ci_lower=stats::quantile(surv_b, probs=1-conf_level,
                                                na.rm=T),
                       ci_upper=stats::quantile(surv_b, probs=conf_level,
                                                na.rm=T),
                       n_boot=sum(!is.na(surv_b)),
                       .groups="drop_last")
  }

  # core of the function
  args <- list(data=data, variable=variable, ev_time=ev_time,
               event=event, conf_int=conf_int, conf_level=conf_level,
               times=times, ...)
  plotdata <- R.utils::doCall(surv_fun, args=args)

  out <- list(adjsurv=plotdata,
              method=method,
              categorical=ifelse(length(unique(data[,variable]))>2, T, F),
              call=match.call())

  if (bootstrap) {
    out$boot_data <- boot_data
    out$boot_data_same_t <- boot_data_same_t
    out$boot_adjsurv <- as.data.frame(boot_stats)
  }

  class(out) <- "adjustedsurv"
  return(out)
}

## perform one bootstrap iteration
adjustedsurv_boot <- function(data, variable, ev_time, event, method,
                              times_input, times, i, surv_fun, levs, ...) {

  indices <- sample(x=rownames(data), size=nrow(data), replace=T)
  boot_samp <- data[indices,]
  # IMPORTANT: keeps SL in tmle methods from failing
  row.names(boot_samp) <- 1:nrow(data)

  # if event specific times are used, use event specific times
  # in bootstrapping as well
  if (is.null(times_input)) {
    times_boot <- sort(unique(boot_samp[, ev_time][boot_samp[, event]==1]))

    if (!0 %in% times_boot) {
      times_boot <- c(0, times_boot)
    }
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
    if (inherits(pass_args$treatment_model, "glm") |
        inherits(pass_args$treatment_model, "multinom")) {
      pass_args$treatment_model <- stats::update(pass_args$treatment_model,
                                                 data=boot_samp, trace=F)
    }
  }

  # call surv_method with correct arguments
  args <- list(data=boot_samp, variable=variable, ev_time=ev_time,
               event=event, conf_int=F, conf_level=0.95, times=times)
  args <- c(args, pass_args)

  adjsurv_boot <- R.utils::doCall(surv_fun, args=args)
  adjsurv_boot$boot <- i

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
    boot_surv[[j]] <- dat_temp
  }
  boot_surv <- as.data.frame(dplyr::bind_rows(boot_surv))

  # output
  result <- multi_result_class()

  result$boot_data <- adjsurv_boot
  result$boot_data_same_t <- boot_surv

  return(result)

}

## plot the survival curves
#' @importFrom rlang .data
#' @export
plot.adjustedsurv <- function(x, draw_ci=F, max_t=Inf,
                              iso_reg=F, force_bounds=F, use_boot=F,
                              color=T, linetype=F, facet=F,
                              line_size=1, xlab="Time",
                              ylab="Adjusted Survival Probability",
                              title=NULL, legend_title="Group",
                              legend_position="right",
                              ylim=NULL, custom_colors=NULL,
                              custom_linetypes=NULL,
                              ci_draw_alpha=0.4, ...) {

  if (!color & !linetype & !facet) {
    stop("Groups must be distinguished with at least one of 'color',",
         "'linetype' or 'facet'. Can't all be FALSE.")
  }

  if (use_boot & is.null(x$boot_adjsurv)) {
    warning("Cannot use bootstrapped estimates as they were not estimated.",
            " Need bootstrap=TRUE in adjustedsurv() call.")
    draw_ci <- F
    plotdata <- x$adjsurv
  } else if (use_boot) {
    plotdata <- x$boot_adjsurv
  } else {
    plotdata <- x$adjsurv
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

  if (draw_ci & !"ci_lower" %in% colnames(plotdata)) {
    warning("Cannot draw confidence intervals. Either set 'conf_int=TRUE' in",
            " 'adjustedsurv()' call or use bootstrap estimates.")
  }

  if (draw_ci & "ci_lower" %in% colnames(plotdata)) {
    p <- p + pammtools::geom_stepribbon(ggplot2::aes(ymin=.data$ci_lower,
                                            ymax=.data$ci_upper,
                                            fill=.data$group,
                                            x=.data$time,
                                            y=.data$surv),
                                        alpha=ci_draw_alpha, inherit.aes=F)
  }
  return(p)
}

## function to calculate the restricted mean survival time of each
## adjusted survival curve previously estimated using the adjustedsurv function
#' @export
adjusted_rmst <- function(adjsurv, to, from=0, use_boot=F, conf_level=0.95) {

  check_inputs_adj_rmst(adjsurv=adjsurv, from=from, to=to, use_boot=use_boot)

  if (use_boot) {

    n_boot <- max(adjsurv$boot_data$boot)
    booted_rmsts <- vector(mode="list", length=n_boot)

    for (i in 1:n_boot) {

      # select one bootstrap data set each
      boot_dat <- adjsurv$boot_data[adjsurv$boot_data$boot==i,]

      # create fake adjustedsurv object
      fake_adjsurv <- list(adjsurv=boot_dat)
      class(fake_adjsurv) <- "adjustedsurv"

      # recursion call
      adj_rmst <- adjusted_rmst(fake_adjsurv, from=from, to=to, use_boot=F)

      booted_rmsts[[i]] <- adj_rmst$rmsts
    }
    booted_rmsts <- as.data.frame(dplyr::bind_rows(booted_rmsts))
  }

  levs <- unique(adjsurv$adjsurv$group)

  rmsts <- vector(mode="numeric", length=length(levs))
  for (i in 1:length(levs)) {

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

  if (use_boot) {

    n_boot_rmst <- apply(booted_rmsts, 2, function(x){sum(!is.na(x))})
    names(n_boot_rmst) <- levs

    out$conf_level <- conf_level
    out$n_boot <- n_boot_rmst
    out$booted_rmsts <- booted_rmsts
    out$rmsts_sd <- apply(booted_rmsts, 2, stats::sd, na.rm=T)
    out$rmsts_ci_lower <- apply(booted_rmsts, 2, stats::quantile,
                                probs=1-conf_level, na.rm=T)
    out$rmsts_ci_upper <- apply(booted_rmsts, 2, stats::quantile,
                                probs=conf_level, na.rm=T)
  }

  return(out)
}

## print method for adjusted_rmst function
#' @export
print.adjusted_rmst <- function(x, digits=5, ...) {

  cat("------------------------------------------------------------------\n")
  cat("Confounder-Adjusted Restricted Mean Survival Time\n")
  cat("------------------------------------------------------------------\n")
  cat("\n")
  cat("Using the interval:", x$from, "to", x$to, "\n")
  cat("\n")

  if (!is.null(x$booted_rmsts)) {
    all_data <- rbind(x$rmsts, x$rmsts_sd, x$rmsts_ci_lower,
                      x$rmsts_ci_upper, x$n_boot)

    ci_lower_name <- paste0(round(x$conf_level*100, 2), "% CI (lower)")
    ci_upper_name <- paste0(round(x$conf_level*100, 2), "% CI (upper)")

    rownames(all_data) <- c("RMST", "RMST SD", ci_lower_name, ci_upper_name,
                            "N Boot")
  } else {
    all_data <- data.frame(x$rmsts)
    rownames(all_data) <- c("RMST")
  }

  colnames(all_data) <- paste0("Group=", colnames(all_data))
  all_data <- round(all_data, digits)
  print(t(all_data))

  cat("------------------------------------------------------------------\n")

  # also silently return that data.frame
  return(invisible(t(all_data)))
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
  group_p <- intercept + rowSums(treatment_betas *
                                 covars[,names(treatment_betas)])
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
