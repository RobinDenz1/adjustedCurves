
# TODO:
#   - input checks for adjust_vars, instrument and .RES.
#   - unit tests
#   - update every man page with estimand information

## estimate counterfactual survival curves using a two-stage instrumental
## variable approach as described in Martinez-Camblor (2021)
#' @export
surv_iv_2SRIF <- function(data, variable, ev_time, event, conf_int=FALSE,
                          conf_level=0.95, times, adjust_vars, instrument,
                          frailty_dist="gaussian", return_models=TRUE) {

  # coerce variable to numeric
  levs <- levels(data[, variable])
  data[, variable] <- ifelse(data[, variable]==levs[1], 0, 1)

  # fit first stage linear model
  res_form <- paste0(variable, " ~ ", instrument, " + ",
                     paste0(adjust_vars, collapse=" + "))
  lm_mod <- stats::lm(formula=stats::as.formula(res_form),
                      data=data)
  data$.RES. <- lm_mod$residuals

  # fit 2SRI-F Cox model
  n <- nrow(data)
  cox_form <- paste0("Surv(", ev_time, ", ", event, ") ~ ",
                     variable, " + .RES. + ",
                     paste0(adjust_vars, collapse=" + "),
                     " + frailty(1:n, distribution='", frailty_dist, "')")
  cox_tsrif <- survival::coxph(formula=stats::as.formula(cox_form), data=data)

  # get predictor under both variable = 0 and variable = 1
  rel_cols <- as.matrix(data[, c(".RES.", adjust_vars)])

  coefs <- cox_tsrif$coefficients
  coefs <- coefs[2:length(coefs)]

  pred_0 <- rowSums(rel_cols %*% diag(coefs)) + cox_tsrif$frail
  pred_1 <- cox_tsrif$coef[1] + pred_0

  # get hazards
  b0 <- survival::basehaz(cox_tsrif, centered=FALSE)
  haz <- stats::approxfun(c(0, b0$time), c(0, b0$hazard))(times)

  f0 <- 0
  f1 <- 0
  for (i in seq_len(n)) {
    f0 <- f0 + exp(-b0$hazard * exp(pred_0[i]))
    f1 <- f1 + exp(-b0$hazard * exp(pred_1[i]))
  }

  # estimate final survival probabilties
  surv_0 <- stats::approxfun(c(0, b0$time), c(1,f0/n))(times)
  surv_1 <- stats::approxfun(c(0, b0$time), c(1,f1/n))(times)

  # put together
  plotdata <- data.frame(time=c(times, times),
                         group=c(rep(0, length(times)),
                                 rep(1, length(times))),
                         surv=c(surv_0, surv_1))

  # coerce variable back to factor
  plotdata$group <- factor(ifelse(plotdata$group==0, levs[1], levs[2]),
                           levels=levs)

  out <- list(plotdata=plotdata)
  class(out) <- "adjustedsurv.method"

  if (return_models) {
    out$lm_mod <- lm_mod
    out$cox_mod <- cox_tsrif
  }

  return(out)
}
