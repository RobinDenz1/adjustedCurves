
## applies the estimator introduced by Cheng et al. (2022), which combines
## weighting based adjustment for confounding with inverse probability of
## censoring weighting
#' @export
surv_iptw_cens <- function(data, variable, ev_time, event, conf_int=FALSE,
                           conf_level=0.95, times=NULL, treatment_model,
                           weight_method="ps", stabilize=FALSE,
                           trim=FALSE, trim_quantiles=FALSE, censoring_model,
                           coxph_control=list(), ties="efron", ...) {
  # preliminaries
  treatment <- data[[variable]]
  treatment_levels <- levels(treatment)
  n_levels <- length(treatment_levels)

  # get weights
  if (is.numeric(treatment_model)) {
    weights <- treatment_model
    weights <- trim_weights(weights=weights, trim=trim)
    weights <- trim_weights_quantiles(weights=weights, trim_q=trim_quantiles)
    if (stabilize) {
      weights <- stabilize_weights(weights, data, variable, treatment_levels)
    }
  } else {
    weights <- get_iptw_weights(data=data, treatment_model=treatment_model,
                                weight_method=weight_method,
                                variable=variable, stabilize=stabilize,
                                trim=trim, trim_q=trim_quantiles, ...)
  }

  # get censoring probabilities
  censoring_result <- fit_censoring_model(
    data=data, time_var=ev_time, event_var=event, treatment_var=variable,
    censoring_formula=censoring_model, control=coxph_control,
    ties=ties
  )

  time <- data[[ev_time]]
  event <- data[[event]]

  # times to evaluate at
  if (is.null(times)) {
    eval_times <- sort(unique(time[event == 1]))
  } else {
    eval_times <- times
  }
  n_times <- length(eval_times)

  survival_matrix <- matrix(
    NA_real_,
    nrow = n_times,
    ncol = n_levels,
    dimnames = list(NULL, treatment_levels)
  )

  weighted_n_risk <- matrix(
    NA_real_,
    nrow = n_times,
    ncol = n_levels,
    dimnames = list(NULL, treatment_levels)
  )

  n_risk <- matrix(
    NA_real_,
    nrow = n_times,
    ncol = n_levels,
    dimnames = list(NULL, treatment_levels)
  )

  weighted_n_events <- matrix(
    0,
    nrow = n_times,
    ncol = n_levels,
    dimnames = list(NULL, treatment_levels)
  )

  n_events <- matrix(
    0,
    nrow = n_times,
    ncol = n_levels,
    dimnames = list(NULL, treatment_levels)
  )

  weighted_n_events_cum <- matrix(
    0,
    nrow = n_times,
    ncol = n_levels,
    dimnames = list(NULL, treatment_levels)
  )

  n_censored <- matrix(
    0,
    nrow = n_times,
    ncol = n_levels,
    dimnames = list(NULL, treatment_levels)
  )

  weighted_n_censored <- matrix(
    0,
    nrow = n_times,
    ncol = n_levels,
    dimnames = list(NULL, treatment_levels)
  )

  for (j in seq_along(treatment_levels)) {

    level <- treatment_levels[j]
    idx <- treatment == level

    w <- weights[idx]
    t <- time[idx]
    d <- event[idx]
    G <- censoring_result$censoring_prob[idx, j]

    denominator <- sum(w)

    # IPCW contribution to the estimated CDF
    event_weight <- w * d / G

    # event times and their weighted contributions
    event_times <- t[d == 1]
    event_weights <- event_weight[d == 1]

    for (k in seq_along(eval_times)) {

      tau <- eval_times[k]

      # Cheng estimator I
      weighted_n_events_cum[k, j] <- sum(event_weight[t <= tau])
      survival_matrix[k, j] <- 1 - weighted_n_events_cum[k, j] / denominator

      # ordinary treatment-weighted risk set
      at_risk <- t >= tau
      n_risk[k, j] <- sum(at_risk)
      weighted_n_risk[k, j] <- sum(w[at_risk])

      # events exactly at this time
      event_at_tau <- d == 1 & t == tau
      n_events[k, j] <- sum(event_at_tau)
      weighted_n_events[k, j] <- sum(event_weight[event_at_tau])

      # censoring exactly at this time
      censored_at_tau <- d == 0 & t == tau
      n_censored[k, j] <- sum(censored_at_tau)
      weighted_n_censored[k, j] <- sum(w[censored_at_tau])
    }
  }

  # put together
  plotdata <- data.frame(time=rep(eval_times, n_levels),
                         surv=c(survival_matrix),
                         group=rep(treatment_levels, each=length(eval_times)))

  n_at_risk <- data.frame(time=rep(eval_times, n_levels),
                          group=rep(treatment_levels, each=length(eval_times)),
                          n_at_risk=c(weighted_n_risk),
                          n_events=c(weighted_n_events))

  output <- list(plotdata=plotdata,
                 weights=weights,
                 n_at_risk=n_at_risk,
                 censoring_models=censoring_result$models)
  class(output) <- "adjustedsurv.method"

  return(output)
}

## fit Cox based censoring models any extract required time-varying
## censoring probabilities for IPCW adjustment
fit_censoring_model <- function(data, time_var, event_var, treatment_var,
                                censoring_formula, control=list(),
                                ties="efron") {

  treatment <- data[[treatment_var]]
  treatment_levels <- levels(treatment)
  n <- nrow(data)

  # create actual censoring formula
  rhs <- paste(deparse(censoring_formula[[2L]]), collapse = "")
  formula <- stats::as.formula(paste0("survival::Surv(", time_var, ", 1 - ",
                                      event_var, ") ~ ", rhs))

  # initialize output objects
  models <- vector("list", length(treatment_levels))
  baseline_hazards <- vector("list", length(treatment_levels))
  linear_predictors <- matrix(NA_real_, nrow=n, ncol=length(treatment_levels),
                              dimnames=list(NULL, treatment_levels))
  censoring_prob <- linear_predictors

  for (j in seq_along(treatment_levels)) {

    level <- treatment_levels[j]
    idx <- treatment == level
    data_j <- data[idx, , drop = FALSE]

    # fit Cox model for censoring for treatment group j
    fit <- survival::coxph(formula=formula, data=data_j, ties=ties, model=TRUE,
                           control=control)

    # estimate non-parametric baseline hazard
    bh <- survival::basehaz(fit, centered=FALSE)

    # estimate censoring probabilities for all required t
    lp <- stats::predict(fit, newdata=data, type="lp", reference="zero")

    models[[j]] <- fit
    baseline_hazards[[j]] <- bh
    linear_predictors[, j] <- lp

    # H_0(T_i)
    k <- findInterval(data[[time_var]], bh$time)

    H0 <- numeric(n)
    keep <- k > 0
    H0[keep] <- bh$hazard[k[keep]]

    # G(T_i | A = j, X_i)
    censoring_prob[, j] <- exp(-H0 * exp(lp))
  }

  names(models) <- treatment_levels
  names(baseline_hazards) <- treatment_levels

  out <- list(
    models = models,
    baseline_hazards = baseline_hazards,
    linear_predictors = linear_predictors,
    censoring_prob = censoring_prob
  )

  return(out)
}
