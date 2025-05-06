
## Hypothesis-Test for the difference between two adjusted survival curves
## or two adjusted cumulative incidence functions
#' @export
adjusted_curve_test <- function(adj, to, from=0, conf_level=0.95,
                                interpolation="steps",
                                group_1=NULL, group_2=NULL) {

  # silence devtools::check() notes
  . <- comparison <- area_est <- p_val <- n_boot <- NULL

  est <- ifelse(inherits(adj, "adjustedsurv"), "surv", "cif")
  adj_method <- adj$method

  treat_labs <- levels(adj$adj$group)

  ## using multiply imputed results
  if (!is.null(adj$mids_analyses)) {

    output <- adjusted_curve_test.MI(adj=adj, to=to, from=from,
                                     conf_level=conf_level, est=est,
                                     interpolation=interpolation,
                                     adj_method=adj_method,
                                     treat_labs=treat_labs,
                                     group_1=group_1, group_2=group_2)
  ## using regular results
  } else {

    check_inputs_adj_test(adj=adj, from=from, to=to)

    # using just two groups in custom order
    if (!is.null(group_1) & !is.null(group_2)) {

      # keep only the two groups, change their level order and call
      # usual adjusted_curve_test function
      adj$adj <- adj$adj[adj$adj$group %in% c(group_1, group_2), ]
      adj$adj$group <- factor(adj$adj$group, levels=c(group_1, group_2))

      adj$boot_data <- adj$boot_data[adj$boot_data$group %in%
                                     c(group_1, group_2), ]
      adj$boot_data$group <- factor(adj$boot_data$group,
                                    levels=c(group_1, group_2))
      adj$categorical <- FALSE

      out <- adjusted_curve_test(adj=adj, to=to, from=from,
                                 conf_level=conf_level,
                                 interpolation=interpolation,
                                 group_1=NULL, group_2=NULL)

    # just two treatments, standard procedure
    } else if (!adj$categorical) {

      # calculate the integral of the difference for every bootstrap sample
      stats_vec <- vector(mode="numeric", length=max(adj$boot_data$boot))
      curve_list <- vector(mode="list", length=max(adj$boot_data$boot))

      for (i in seq_len(max(adj$boot_data$boot))) {

        # select one bootstrap data set each
        boot_dat <- adj$boot_data[adj$boot_data$boot==i, ]

        # 1.) get every relevant point in time
        # 2.) create new curve of the difference
        # 3.) get integral of that curve
        times <- sort(unique(boot_dat$time))
        surv_diff <- difference_function(adj=boot_dat, times=times, est=est,
                                         interpolation=interpolation,
                                         conf_int=FALSE, conf_level=0.95)
        diff_integral <- exact_integral(data=surv_diff, from=from, to=to,
                                        est=est, interpolation=interpolation)

        stats_vec[i] <- diff_integral

        # also append difference curve to list
        surv_diff$boot <- i
        curve_list[[i]] <- surv_diff
      }

      diff_curves <- as.data.frame(dplyr::bind_rows(curve_list))

      ## use those bootstrapped integrals for the calculation of
      ## a p-value, by shifting the observed distribution to 0 and
      ## comparing it to the actually observed value

      # actually observed values
      adj_observed <- adj$adj

      observed_diff_curve <- difference_function(adj=adj_observed, times=times,
                                                 est=est,
                                                 interpolation=interpolation,
                                                 conf_int=FALSE,
                                                 conf_level=0.95)
      observed_diff_integral <- exact_integral(data=observed_diff_curve,
                                               from=from, to=to, est=est,
                                               interpolation=interpolation)

      # remove NA values
      stats_vec <- stats_vec[!is.na(stats_vec)]

      # shift bootstrap distribution
      diff_under_H0 <- stats_vec - mean(stats_vec)
      p_value <- mean(abs(diff_under_H0) > abs(observed_diff_integral))

      # bootstrap confidence interval
      conf_int <- stats::quantile(stats_vec,
                                  probs=c((1-conf_level)/2,
                                          1-((1-conf_level)/2)),
                                  na.rm=TRUE)
      names(conf_int) <- c("ci_lower", "ci_upper")

      # do this so it shows up in print function
      fun_call <- match.call()
      fun_call$from <- from
      fun_call$to <- to

      # put together output
      out <- list(diff_curves=diff_curves,
                  diff_integrals=stats_vec,
                  observed_diff_curve=observed_diff_curve,
                  observed_diff_integral=observed_diff_integral,
                  integral_se=stats::sd(stats_vec),
                  p_value=p_value,
                  n_boot=length(stats_vec),
                  kind=est,
                  conf_int=conf_int,
                  categorical=FALSE,
                  treat_labs=treat_labs,
                  method=adj_method,
                  interpolation=interpolation,
                  call=fun_call)

    ## more than two treatments -> perform pairwise comparisons
    } else {

      if (inherits(adj, "adjustedsurv")) {
        combs <- all_combs_length_2(treat_labs)
      } else {
        combs <- all_combs_length_2(treat_labs)
      }

      out <- list()
      for (i in seq_len(length(combs))) {

        # get first and second group
        group_0 <- strsplit(combs[i], "\t")[[1]][1]
        group_1 <- strsplit(combs[i], "\t")[[1]][2]

        # create pseudo adjustedsurv, adjustedcif object
        observed_dat <- adj$adj[which(adj$adj$group %in%
                                        c(group_0, group_1)), ]
        observed_dat$group <- factor(observed_dat$group,
                                     levels=c(group_0, group_1))

        boot_dat <- adj$boot_data[which(adj$boot_data$group %in%
                                          c(group_0, group_1)), ]
        boot_dat$group <- factor(boot_dat$group,
                                 levels=c(group_0, group_1))

        fake_adjsurv <- list(adj=observed_dat,
                             boot_data=boot_dat,
                             categorical=FALSE,
                             method=adj_method)
        class(fake_adjsurv) <- class(adj)

        # recursion call
        pair <- adjusted_curve_test(adj=fake_adjsurv,
                                    from=from,
                                    to=to,
                                    conf_level=conf_level,
                                    interpolation=interpolation,
                                    group_1=group_1, group_2=group_2)
        out[[paste0(group_0, " vs. ", group_1)]] <- pair
      }
      out$categorical <- TRUE
      out$method <- adj_method
    }
    class(out) <- "curve_test"

    return(out)
  }
}

## Same test but using multiple imputation
adjusted_curve_test.MI <- function(adj, to, from, conf_level, est,
                                   adj_method, treat_labs,
                                   interpolation, group_1, group_2) {

  # silence devtools::check() notes
  . <- comparison <- area_est <- p_val <- n_boot <- NULL

  if (adj$categorical) {

    # call function once on every adjsurv object, extract values
    len <- length(adj$mids_analyses)
    mids_out <- dat <- vector(mode="list", length=len)
    for (i in seq_len(len)) {

      results_imp <- adjusted_curve_test(adj$mids_analyses[[i]],
                                         to=to, from=from,
                                         conf_level=conf_level,
                                         interpolation=interpolation,
                                         group_1=group_1, group_2=group_2)
      mids_out[[i]] <- results_imp
      comp_names <- names(results_imp)

      # for each pairwise comparison, extract values
      for (j in seq_len((length(results_imp)-2))) {

        row <- data.frame(area_est=results_imp[[j]]$observed_diff_integral,
                          area_se=results_imp[[j]]$integral_se,
                          p_val=results_imp[[j]]$p_value,
                          n_boot=results_imp[[j]]$n_boot,
                          comparison=comp_names[j])
        dat[[length(dat)+1]] <- row
      }

    }
    dat <- dplyr::bind_rows(dat)

    # pool values
    pooled_dat <- dat %>%
      dplyr::group_by(., comparison) %>%
      dplyr::summarise(area_est=mean(area_est),
                       area_se=mean(area_se),
                       p_value=pool_p_values(p_val),
                       n_boot=min(n_boot))

    # create one curve_test object for each comparison
    output <- list(mids_analyses=mids_out)
    for (i in seq_len(nrow(pooled_dat))) {

      pooled_ci <- confint_surv(surv=pooled_dat$area_est[i],
                                se=pooled_dat$area_se[i],
                                conf_level=conf_level,
                                conf_type="plain")
      pooled_ci <- c(pooled_ci$left, pooled_ci$right)
      names(pooled_ci) <- c("ci_lower", "ci_upper")

      # do this so it shows up in print function
      fun_call <- match.call()
      fun_call$from <- from
      fun_call$to <- to

      mids_p <- dat$p_val[pooled_dat$comparison[i]==dat$comparison]
      out <- list(mids_p_values=mids_p,
                  observed_diff_integral=pooled_dat$area_est[i],
                  p_value=pooled_dat$p_value[i],
                  n_boot=pooled_dat$n_boot[i],
                  kind=est,
                  integral_se=pooled_dat$area_se[i],
                  conf_int=pooled_ci,
                  categorical=FALSE,
                  treat_labs=mids_out[[1]][[i]]$treat_labs,
                  method=adj_method,
                  interpolation=interpolation,
                  call=fun_call)
      class(out) <- "curve_test"
      output[[pooled_dat$comparison[i]]] <- out

    }
    output$categorical <- TRUE
    output$method <- adj_method
    class(output) <- "curve_test"

    return(output)

  } else {

    len <- length(adj$mids_analyses)
    mids_out <- vector(mode="list", length=len)
    area_ests <- area_se <- p_vals <- n_boots <- vector(mode="numeric",
                                                        length=len)
    for (i in seq_len(len)) {

      results_imp <- adjusted_curve_test(adj$mids_analyses[[i]],
                                         to=to, from=from,
                                         conf_level=conf_level,
                                         interpolation=interpolation,
                                         group_1=group_1, group_2=group_2)
      mids_out[[i]] <- results_imp
      area_ests[i] <- results_imp$observed_diff_integral
      area_se[i] <- results_imp$integral_se
      p_vals[i] <- results_imp$p_value
      n_boots[i] <- results_imp$n_boot

    }

    # pool the values
    observed_diff_integral <- mean(area_ests)
    pooled_se <- mean(area_se)
    pooled_ci <- confint_surv(surv=observed_diff_integral,
                              se=pooled_se,
                              conf_level=conf_level,
                              conf_type="plain")
    pooled_ci <- c(pooled_ci$left, pooled_ci$right)
    names(pooled_ci) <- c("ci_lower", "ci_upper")

    pooled_p_value <- pool_p_values(p_values=p_vals)

    # do this so it shows up in print function
    fun_call <- match.call()
    fun_call$from <- from
    fun_call$to <- to

    out <- list(mids_analyses=mids_out,
                mids_p_values=p_vals,
                observed_diff_integral=observed_diff_integral,
                p_value=pooled_p_value,
                n_boot=min(n_boots),
                kind=est,
                integral_se=pooled_se,
                conf_int=pooled_ci,
                categorical=FALSE,
                treat_labs=treat_labs,
                method=adj_method,
                interpolation=interpolation,
                call=fun_call)
    class(out) <- "curve_test"

  }
  return(out)
}

## create data.frame of pairwise comparisons
gather_pairwise_comps <- function(x, conf_level) {

  if (!is.null(x$mids_analyses)) {
    sequence <- seq(2, (length(x)-2))
  } else {
    sequence <- seq(1, (length(x)-2))
  }

  x_names <- names(x)
  out <- vector(mode="list", length=(length(x)-2))
  for (i in sequence) {

    out[[i]] <- data.frame(comparison=x_names[[i]],
                           ABC=x[[i]]$observed_diff_integral,
                           ABC_SE=x[[i]]$integral_se,
                           ci_lower=x[[i]]$conf_int[1],
                           ci_upper=x[[i]]$conf_int[2],
                           p_value=x[[i]]$p_value,
                           n_boot=x[[i]]$n_boot)
  }
  out <- as.data.frame(dplyr::bind_rows(out))
  rownames(out) <- out$comparison
  out$comparison <- NULL
  colnames(out) <- c("ABC", "ABC SE",
                     paste0(conf_level, "% CI (lower)"),
                     paste0(conf_level, "% CI (upper)"),
                     "P-Value", "N Boot")

  return(out)
}

## print method for curve_test objects
#' @export
print.curve_test <- function(x, digits=4, ...) {

  if (is.null(x$kind)) {
    kind <- x[[2]]$kind
  } else {
    kind <- x$kind
  }

  if (kind=="surv" & x$method=="km") {
    title <- "   Test of the Difference between two Survival Curves\n"
  } else if (kind=="surv") {
    title <- "   Test of the Difference between two adjusted Survival Curves\n"
  } else if (kind=="cif" & x$method=="aalen_johansen") {
    title <- "   Test of the Difference between two CIFs \n"
  } else if (kind=="cif") {
    title <- "   Test of the Difference between two adjusted CIFs \n"
  }

  if (is.null(x$call)) {
    call_conf <- x[[2]]$call$conf_level
  } else {
    call_conf <- x$call$conf_level
  }

  if (is.null(call_conf) || call_conf=="conf_level") {
    conf_level <- 95
  } else {
    conf_level <- call_conf * 100
  }

  if (x$categorical) {

    out <- gather_pairwise_comps(x=x, conf_level=conf_level)

    cat("------------------------------------------------------------------\n")
    cat(title)
    cat("------------------------------------------------------------------\n")
    cat("\n")
    cat("Using the interval:", x[[2]]$call$from, "to", x[[2]]$call$to, "\n")
    cat("\n")
    print(round(out, digits), row.names=TRUE)
    cat("------------------------------------------------------------------\n")

  } else {

    out <- data.frame(ABC=x$observed_diff_integral,
                      ABC_SE=x$integral_se,
                      ci_lower=x$conf_int[1],
                      ci_upper=x$conf_int[2],
                      p_value=x$p_value,
                      n_boot=x$n_boot)
    colnames(out) <- c("ABC", "ABC SE",
                       paste0(conf_level, "% CI (lower)"),
                       paste0(conf_level, "% CI (upper)"),
                       "P-Value", "N Boot")
    rownames(out) <- paste0(x$treat_labs[1], " vs. ", x$treat_labs[2])

    cat("------------------------------------------------------------------\n")
    cat(title)
    cat("------------------------------------------------------------------\n")
    cat("\n")
    cat("Using the interval:", x$call$from, "to", x$call$to, "\n")
    cat("\n")
    print(round(out, digits), row.names=TRUE)
    cat("------------------------------------------------------------------\n")
  }

  # also silently return the data.frame
  return(invisible(out))
}

## summary method for curve_test objects
#' @export
summary.curve_test <- function(object, ...) {
  print.curve_test(object, ...)
}

## plot method for curve_test objects
#' @export
plot.curve_test <- function(x, type="curves", xlab=NULL, ylab=NULL,
                            title=NULL, ...) {
  requireNamespace("ggplot2")

  # to remove devtools::check() Notes
  time <- surv <- boot <- integral <- NULL

  if (!is.null(x$mids_analyses)) {

    stop("There is no plot method for 'curve_test' objects fitted using",
         " multiply imputed datasets.")

  } else if (x$categorical) {

    comps <- names(x)

    if (type=="integral") {

      observed_diff_integrals <- list()
      diff_integrals <- list()
      for (i in seq_len((length(x)-2))) {

        # observed diff curves
        obs_diff <- x[[i]]$observed_diff_integral
        observed_diff_integrals[[i]] <- data.frame(integral=obs_diff,
                                                   comp=comps[i])

        # bootstrapped diff curves
        diff <- x[[i]]$diff_integrals

        # remove NA values
        diff <- diff[!is.na(diff)]

        # shift bootstrap distribution
        diff <- diff - mean(diff)

        diff_integrals[[i]] <- data.frame(integral=diff,
                                          comp=comps[i])
      }
      observed_diff_integrals <- dplyr::bind_rows(observed_diff_integrals)
      diff_integrals <- dplyr::bind_rows(diff_integrals)

      if (is.null(xlab)) {
        xlab <- "Integrals of the difference under H0"
      }

      if (is.null(ylab)) {
        ylab <- "Density"
      }

      # plot it
      p <- ggplot2::ggplot(diff_integrals, ggplot2::aes(x=integral)) +
        ggplot2::geom_density() +
        ggplot2::geom_vline(xintercept=0, linetype="dashed") +
        ggplot2::geom_vline(data=observed_diff_integrals,
                            ggplot2::aes(xintercept=integral),
                            color="red") +
        ggplot2::theme_bw() +
        ggplot2::labs(x=xlab, y=ylab) +
        ggplot2::facet_wrap(~comp, scales="free")

    } else if (type=="curves") {

      if (is.null(ylab)) {
        ylab <- ifelse(x[[1]]$kind=="surv", "Difference in Survival",
                       "Difference in Cumulative Incidence")
      }

      if (is.null(xlab)) {
        xlab <- "Time"
      }

      observed_diff_curves <- list()
      diff_curves <- list()
      for (i in seq_len((length(x)-2))) {

        # observed diff curves
        obs_diff <- x[[i]]$observed_diff_curve
        obs_diff$comp <- comps[i]
        observed_diff_curves[[i]] <- obs_diff

        # bootstrapped diff curves
        diff <- x[[i]]$diff_curves
        diff$comp <- comps[i]
        diff_curves[[i]] <- diff
      }
      observed_diff_curves <- dplyr::bind_rows(observed_diff_curves)
      colnames(observed_diff_curves) <- c("time", "surv", "comp")
      diff_curves <- dplyr::bind_rows(diff_curves)
      colnames(diff_curves) <- c("time", "surv", "boot", "comp")

      p <- ggplot2::ggplot(diff_curves, ggplot2::aes(x=time, y=surv))

      if (x[[1]]$interpolation=="steps") {
        p <- p + ggplot2::geom_step(ggplot2::aes(group=boot), color="grey",
                                    alpha=0.8) +
          ggplot2::geom_step(data=observed_diff_curves,
                             ggplot2::aes(x=time, y=surv))
      } else {
        p <- p + ggplot2::geom_line(ggplot2::aes(group=boot), color="grey",
                                    alpha=0.8) +
          ggplot2::geom_line(data=observed_diff_curves,
                             ggplot2::aes(x=time, y=surv))
      }

      p <- p +
        ggplot2::geom_hline(yintercept=0, linetype="dashed") +
        ggplot2::theme_bw() +
        ggplot2::labs(x=xlab, y=ylab) +
        ggplot2::facet_wrap(~comp, scales="free")
    } else {
      stop("type='", type, "' is not defined.")
    }

  } else {

    if (type=="integral") {

      # remove NA values
      stats_vec <- x$diff_integrals[!is.na(x$diff_integrals)]

      # shift bootstrap distribution
      diff_under_H0 <- stats_vec - mean(stats_vec)

      if (is.null(xlab)) {
        xlab <- "Integrals of the difference under H0"
      }

      if (is.null(ylab)) {
        ylab <- "Density"
      }

      p <- ggplot2::ggplot(NULL, ggplot2::aes(x=diff_under_H0)) +
        ggplot2::geom_density() +
        ggplot2::geom_vline(xintercept=0, linetype="dashed") +
        ggplot2::geom_vline(xintercept=x$observed_diff_integral, color="red") +
        ggplot2::theme_bw() +
        ggplot2::labs(x=xlab, y=ylab)

    } else if (type=="curves") {

      if (is.null(ylab)) {
        ylab <- ifelse(x$kind=="surv", "Difference in Survival",
                       "Difference in Cumulative Incidence")
      }

      if (is.null(xlab)) {
        xlab <- "Time"
      }

      diff_curves <- x$diff_curves
      colnames(diff_curves) <- c("time", "surv", "boot")

      observed_diff_curve <- x$observed_diff_curve
      colnames(observed_diff_curve) <- c("time", "surv")

      p <- ggplot2::ggplot(diff_curves, ggplot2::aes(x=time, y=surv))

      if (x$interpolation=="steps") {
        p <- p + ggplot2::geom_step(ggplot2::aes(group=boot), color="grey",
                                    alpha=0.8) +
          ggplot2::geom_step(data=observed_diff_curve,
                             ggplot2::aes(x=time, y=surv))
      } else {
        p <- p + ggplot2::geom_line(ggplot2::aes(group=boot), color="grey",
                                    alpha=0.8) +
          ggplot2::geom_line(data=observed_diff_curve,
                             ggplot2::aes(x=time, y=surv))
      }

      p <- p + ggplot2::geom_hline(yintercept=0, linetype="dashed") +
        ggplot2::theme_bw() +
        ggplot2::labs(x=xlab, y=ylab)

    } else {
      stop("type='", type, "' is not defined.")
    }
  }

  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }

  return(p)
}

## function that calculates all possible combinations
## of length 2 without replacement and without order
all_combs_length_2 <- function(treat_labs) {

  combs <- list()
  for (i in seq_len(length(treat_labs))) {
    for (j in seq_len(length(treat_labs))) {

      # no cases with the same group twice and order does not matter
      if (i != j & !paste(treat_labs[j], treat_labs[i], sep="\t") %in% combs) {
        combs[[length(combs)+1]] <- paste(treat_labs[i],
                                          treat_labs[j], sep="\t")
      }
    }
  }
  combs <- unlist(combs)
  return(combs)
}

## directly pool the p-values using the method by Licht & Rubin
pool_p_values <- function(p_values, tol=1e-14) {

  # if they are all equal, just return that value
  if (stats::var(p_values)==0) {
    return(p_values[1])
  }

  # apply tolerance to p values that are 0
  p_values[p_values==0] <- p_values[p_values==0] + tol

  # transform to z-scale
  z <- stats::qnorm(p_values)
  num <- mean(z)
  den <- sqrt(1 + stats::var(z))
  pooled <- stats::pnorm(num / den)

  return(pooled)
}
