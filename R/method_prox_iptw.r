
## utility function to get approximate standard errors for survival
## probability estimates from PIPTW
get_piptw_se <- function(data, variable, ev_time, sorted_time,
                         noncens, q_bridge, surv, surv_sub, a) {

  # get influence function
  surv_IF <- surv_sub - surv / nrow(data)

  q_dot <- PIPW_q_dot(sorted_time=sorted_time,
                      time=data[, ev_time],
                      noncensor_cumhaz=noncens$noncensor_cumhaz,
                      A=data[, variable],
                      q_Z=q_bridge$q_Z,
                      q_bridge$t,
                      a=a)
  q_censor <- PIPW_censor_dot(sorted_time=sorted_time,
                              time=data[, ev_time],
                              noncensor_cumhaz=noncens$noncensor_cumhaz,
                              A=data[, variable],
                              q_Z=q_bridge$q_Z,
                              q_bridge$t,
                              a=a)

  K <- length(sorted_time)
  IF <- matrix(0, nrow=K, ncol=nrow(data))
  for (i in seq_len(nrow(data))) {
    for (j in seq_len(K)) {
      IF[j, i] <- surv_IF[j, i] +
        sum(q_bridge$t_IF[, i] * q_dot[j, ]) +
        sum(noncens$noncensor_cumhaz_IF[, i] * q_censor[j, ])
    }
  }

  out <- sqrt(rowSums(IF^2))

  return(out)
}

## proximal inverse probability of treatment weighting
#' @export
surv_prox_iptw <- function(data, variable, ev_time, event, conf_int,
                           conf_level=0.95, times=NULL, adjust_vars,
                           treatment_proxy, outcome_proxy,
                           optim_method="BFGS", optim_control=list(),
                           return_fit=TRUE) {

  # turn variable into numeric
  levs <- levels(data[, variable])
  data[, variable] <- ifelse(data[, variable]==levs[1], 0, 1)

  # (sorted) times at which events occur
  sorted_time <- sort(unique(data[, ev_time]))

  # cumulative hazard of event
  noncens <- get_noncensor_cumhaz(data=data, ev_time=ev_time, event=event,
                                  sorted_time=sorted_time)

  # q confounding bridge
  q_bridge <- get_q_confounding_bridge(data=data,
                                       variable=variable,
                                       adjust_vars=adjust_vars,
                                       treatment_proxy=treatment_proxy,
                                       outcome_proxy=outcome_proxy,
                                       optim_method=optim_method,
                                       optim_control=optim_control)

  # PIPW estimates for variable = 0
  surv_0_sub <- PIPW_surv(sorted_time=sorted_time,
                          time=data[, ev_time],
                          noncensor_cumhaz=noncens$noncensor_cumhaz,
                          A=data[, variable],
                          q_Z=q_bridge$q_Z,
                          t=q_bridge$t,
                          a=0)
  surv_0 <- rowSums(surv_0_sub)

  # PIPW estimates for variable = 1
  surv_1_sub <- PIPW_surv(sorted_time=sorted_time,
                          time=data[, ev_time],
                          noncensor_cumhaz=noncens$noncensor_cumhaz,
                          A=data[, variable],
                          q_Z=q_bridge$q_Z,
                          t=q_bridge$t,
                          a=1)
  surv_1 <- rowSums(surv_1_sub)

  # put it together
  plotdata <- data.frame(time=c(sorted_time, sorted_time),
                         surv=c(surv_0, surv_1),
                         group=c(rep(levs[1], length(sorted_time)),
                                 rep(levs[2], length(sorted_time))))

  # turn variable back to factor
  plotdata$group <- factor(plotdata$group, levels=levs)

  # estimate approximate standard errors + CI if specified
  if (conf_int) {
    surv_0_se <- get_piptw_se(data=data,
                              variable=variable,
                              ev_time=ev_time,
                              sorted_time=sorted_time,
                              noncens=noncens,
                              q_bridge=q_bridge,
                              surv=surv_0,
                              surv_sub=surv_0_sub,
                              a=0)

    surv_1_se <- get_piptw_se(data=data,
                              variable=variable,
                              ev_time=ev_time,
                              sorted_time=sorted_time,
                              noncens=noncens,
                              q_bridge=q_bridge,
                              surv=surv_1,
                              surv_sub=surv_1_sub,
                              a=1)

    plotdata$se <- c(surv_0_se, surv_1_se)

    # approximate confidence intervals
    surv_ci <- confint_surv(surv=plotdata$surv,
                            se=plotdata$se,
                            conf_level=conf_level,
                            conf_type="plain")
    plotdata$ci_lower <- surv_ci$left
    plotdata$ci_upper <- surv_ci$right
  }

  if (!is.null(times)) {
    plotdata <- specific_times(plotdata, times)
  }

  if (return_fit) {
    output <- list(plotdata=plotdata,
                   noncensor_cumhaz=noncens$noncensor_cumhaz,
                   noncensor_cumhaz_IF=noncens$noncensor_cumhaz_IF,
                   q_bridge=q_bridge$qlink)
  } else {
    output <- list(plotdata=plotdata)
  }

  class(output) <- "adjustedsurv.method"

  return(output)
}
