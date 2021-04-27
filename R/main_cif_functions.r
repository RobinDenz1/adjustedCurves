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

## Second main function of the package. Is basically a wrapper around
## all cif_method functions, offering additional high level stuff
#' @importFrom dplyr %>%
#' @importFrom doRNG %dorng%
#' @importFrom foreach %dopar%
#' @export
adjustedcif <- function(data, variable, ev_time, event, cause, method,
                        conf_int=F, conf_level=0.95, times=NULL, bootstrap=F,
                        n_boot=500, n_cores=1, ...) {

  check_inputs_adjustedcif(data=data, variable=variable, ev_time=ev_time,
                           event=event, cause=cause, method=method,
                           conf_int=conf_int, conf_level=conf_level,
                           times=times, bootstrap=bootstrap,
                           n_boot=n_boot, ...)

  # define those to remove Notes in devtools::check()
  . <- i <- time <- group <- cif_b <- NULL

  # get event specific times
  times_input <- times
  if (is.null(times) & method=="aalen_johansen") {
    times <- NULL
  } else if (is.null(times)) {
    times <- sort(unique(data[, ev_time][data[, event]>=1]))

    # add zero if not already in there
    if (!0 %in% times) {
      times <- c(0, times)
    }
  }

  levs <- unique(data[,variable])

  # get relevant cif_method function
  cif_fun <- get(paste0("cif_", method))

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
                       "multi_result_class", "adjustedcif_boot")

      boot_out <- foreach::foreach(i=1:n_boot, .packages=pkgs,
                                   .export=export_objs) %dorng% {

        adjustedcif_boot(data=data, variable=variable, ev_time=ev_time,
                         event=event, method=method, times_input=times_input,
                         times=times, i=i, cif_fun=cif_fun,
                         levs=levs, cause=cause, ...)
      }
      parallel::stopCluster(cl)

    } else {

      boot_out <- vector(mode="list", length=n_boot)
      for (i in 1:n_boot) {
        boot_out[[i]] <- adjustedcif_boot(data=data, variable=variable,
                                          ev_time=ev_time, event=event,
                                          method=method, times_input=times_input,
                                          times=times, i=i, cause=cause,
                                          cif_fun=cif_fun, levs=levs, ...)
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
      dplyr::summarise(cif=mean(cif_b, na.rm=T),
                       sd=stats::sd(cif_b, na.rm=T),
                       ci_lower=stats::quantile(cif_b,
                                                probs=(1-conf_level)/2,
                                                na.rm=T),
                       ci_upper=stats::quantile(cif_b,
                                                probs=1-((1-conf_level)/2),
                                                na.rm=T),
                       n_boot=sum(!is.na(cif_b)),
                       .groups="drop_last")
  }

  # core of the function
  args <- list(data=data, variable=variable, ev_time=ev_time,
               event=event, conf_int=conf_int, conf_level=conf_level,
               times=times, cause=cause, ...)
  plotdata <- R.utils::doCall(cif_fun, args=args)

  out <- list(adjcif=plotdata,
              method=method,
              categorical=ifelse(length(unique(data[,variable]))>2, T, F),
              call=match.call())

  if (bootstrap) {
    out$boot_data <- boot_data
    out$boot_data_same_t <- boot_data_same_t
    out$boot_adjcif <- as.data.frame(boot_stats)
  }

  class(out) <- "adjustedcif"
  return(out)
}

## perform one bootstrap iteration
adjustedcif_boot <- function(data, variable, ev_time, event, cause, method,
                             times_input, times, i, cif_fun, levs, ...) {

  indices <- sample(x=rownames(data), size=nrow(data), replace=T)
  boot_samp <- data[indices,]
  # IMPORTANT: keeps SL in tmle methods from failing
  row.names(boot_samp) <- 1:nrow(data)

  # if event specific times are used, use event specific times
  # in bootstrapping as well
  if (is.null(times_input) & method=="aalen_johansen") {

  } else if (is.null(times_input)) {
    times_boot <- sort(unique(boot_samp[, ev_time][boot_samp[, event]>=1]))

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

  if (method %in% c("iptw", "iptw_pseudo", "aiptw", "aiptw_pseudo")) {
    if (inherits(pass_args$treatment_model, "glm") |
        inherits(pass_args$treatment_model, "multinom")) {
      pass_args$treatment_model <- stats::update(pass_args$treatment_model,
                                                 data=boot_samp, trace=F)
    }
  }

  # call cif_method with correct arguments
  args <- list(data=boot_samp, variable=variable, ev_time=ev_time,
               event=event, conf_int=F, conf_level=0.95, times=times,
               cause=cause)
  args <- c(args, pass_args)

  adjcif_boot <- R.utils::doCall(cif_fun, args=args)
  adjcif_boot$boot <- i

  # read from resulting step function at all t in times
  boot_cif <- vector(mode="list", length=length(levs))
  for (j in 1:length(levs)) {

    if (method=="aalen_johansen" & is.null(times)) {
      times <- unique(data[,ev_time][data[,variable]==levs[j]])
    }

    cif_boot <- sapply(times, read_from_step_function,
                       step_data=adjcif_boot[adjcif_boot$group==levs[j],],
                       est="cif")

    dat_temp <- data.frame(time=times,
                           cif_b=cif_boot,
                           group=levs[j],
                           boot=i)
    boot_cif[[j]] <- dat_temp
  }
  boot_cif <- as.data.frame(dplyr::bind_rows(boot_cif))

  # output
  result <- multi_result_class()

  result$boot_data <- adjcif_boot
  result$boot_data_same_t <- boot_cif

  return(result)

}


## plot the cumulative incidence functions
#' @importFrom rlang .data
#' @export
plot.adjustedcif <- function(x, draw_ci=F, max_t=Inf,
                             iso_reg=F, force_bounds=F, use_boot=F,
                             color=T, linetype=F, facet=F,
                             line_size=1, xlab="Time",
                             ylab="Adjusted Cumulative Incidence",
                             title=NULL, legend_title="Group",
                             legend_position="right",
                             ylim=NULL, custom_colors=NULL,
                             custom_linetypes=NULL,
                             ci_draw_alpha=0.4, ...) {

  if (!color & !linetype & !facet) {
    stop("Groups must be distinguished with at least one of 'color',",
         "'linetype' or 'facet'. Can't all be FALSE.")
  }

  if (use_boot & is.null(x$boot_adjcif)) {
    warning("Cannot use bootstrapped estimates as they were not estimated.",
            " Need bootstrap=TRUE in adjustedcif() call.")
    draw_ci <- F
    plotdata <- x$adjcif
  } else if (use_boot) {
    plotdata <- x$boot_adjcif
  } else {
    plotdata <- x$adjcif
  }
  plotdata$group <- factor(plotdata$group)

  # shortcut to only show curves up to a certain time
  plotdata <- plotdata[which(plotdata$time <= max_t),]

  # in some methods estimates can be outside the 0, 1 bounds,
  # if specified set those to 0 or 1 respectively
  if (force_bounds) {
    plotdata <- within(plotdata, {
      cif <- ifelse(cif < 0, 0, cif);
      cif <- ifelse(cif > 1, 1, cif);
    })
  }

  # apply isotonic regression if specified
  if (iso_reg) {
    for (lev in levels(plotdata$group)) {
      cif <- plotdata$cif[plotdata$group==lev]
      new <- stats::isoreg(cif)$yf
      plotdata$cif[plotdata$group==lev] <- new
    }
  }

  mapping <- ggplot2::aes(x=.data$time, y=.data$cif, color=.data$group,
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
            " 'adjustedcif()' call or use bootstrap estimates.")
  }

  if (draw_ci & "ci_lower" %in% colnames(plotdata)) {
    p <- p + pammtools::geom_stepribbon(ggplot2::aes(ymin=.data$ci_lower,
                                                     ymax=.data$ci_upper,
                                                     fill=.data$group,
                                                     x=.data$time,
                                                     y=.data$cif),
                                        alpha=ci_draw_alpha, inherit.aes=F)
  }
  return(p)
}

## function to simulate confounded competing risks data
#' @export
sim_confounded_crisk <- function(n=500, lcovars=NULL, outcome_betas=NULL,
                                 group_beta=c(1, 0), gamma=c(1.8, 1.8),
                                 lambda=c(2, 2), treatment_betas=NULL,
                                 intercept=-0.5, gtol=0.001,
                                 cens_fun=function(n){stats::rweibull(n, 1, 2)},
                                 cens_args=list(), max_t=1.7) {

  check_inputs_sim_crisk_fun(n=n, lcovars=lcovars, outcome_betas=outcome_betas,
                             gamma=gamma, lambda=lambda,
                             treatment_betas=treatment_betas,
                             group_beta=group_beta, intercept=intercept,
                             gtol=gtol, cens_fun=cens_fun, cens_args=cens_args,
                             max_t=max_t)

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
    # TODO: set better default values
    outcome_betas <- list(c(log(1.8), log(1.8)/3),
                          c(log(1.3), log(1.3)/3),
                          c(0, 0),
                          c(log(1.8), log(1.8)/3),
                          c(log(1.3), log(1.3)/3),
                          c(0, 0))
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
  if (length(treatment_betas)==1) {
    group_p <- intercept + (treatment_betas * covars[,names(treatment_betas)])
  } else {
    group_p <- intercept + rowSums(treatment_betas *
                                     covars[,names(treatment_betas)])
  }
  group_p <- 1/(1 + exp(-group_p))

  # in order to keep the positivity assumption,
  # values of 0 and 1 can not be tolerated
  group_p[group_p < gtol] <- gtol
  group_p[group_p > (1 - gtol)] <- 1 - gtol

  covars$group <- stats::rbinom(n=n, size=1, prob=group_p)

  # add group_beta to outcome_betas
  outcome_betas[[length(outcome_betas) + 1]] <- group_beta


  # generate cause-specific survival times and cause
  surv <- apply(covars, 1, sim_crisk_time, outcome_betas=outcome_betas,
                gamma=gamma, lambda=lambda, max_t=max_t)
  surv <- as.data.frame(t(surv))

  covars$time <- surv$time
  covars$event <- surv$cause

  # introduce random censoring if specified
  if (!is.null(cens_fun)) {
    cens_time <- do.call(cens_fun, args=c(n=n, cens_args))

    covars <- within(covars, {
      event <- ifelse(time < cens_time, event, 0);
      time <- ifelse(time < cens_time, time, cens_time);
    })
  }

  return(covars)
}
