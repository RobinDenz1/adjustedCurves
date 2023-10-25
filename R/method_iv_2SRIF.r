
## if a formula ends with a plus, remove it
## (may happen if adjust_vars is NULL)
remove_plus_at_end <- function(form_str) {
  if (endsWith(form_str, "+ ")) {
    form_str <- substr(form_str, 1, nchar(form_str)-2)
  }
  return(form_str)
}

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
  res_form <- remove_plus_at_end(res_form)

  lm_mod <- stats::lm(formula=stats::as.formula(res_form),
                      data=data)
  data$.RES. <- lm_mod$residuals

  # fit 2SRI-F Cox model
  n <- nrow(data)
  cox_form <- paste0("Surv(", ev_time, ", ", event, ") ~ ",
                     variable, " + .RES. + ",
                     paste0(adjust_vars, collapse=" + "),
                     " + frailty(1:n, distribution='", frailty_dist, "')")
  cox_form <- remove_plus_at_end(cox_form)

  cox_tsrif <- survival::coxph(formula=stats::as.formula(cox_form), data=data)

  # create design matrix for the case of categorical covariates
  x_dat <- as.data.frame(data[, c(".RES.", adjust_vars)])
  if (ncol(x_dat)==1) {
    colnames(x_dat) <- ".RES."
  }
  form <- paste0("~ ", paste0(c(".RES.", adjust_vars), collapse=" + "))
  form <- remove_plus_at_end(form)

  mod_mat <- stats::model.matrix(stats::as.formula(form), data=x_dat)
  mod_mat <- as.matrix(mod_mat[,seq(2, ncol(mod_mat))])

  # get predictor under both variable = 0 and variable = 1
  coefs <- cox_tsrif$coefficients
  coefs <- coefs[2:length(coefs)]

  if (ncol(mod_mat)==1) {
    pred_0 <- mod_mat[,1]*coefs + cox_tsrif$frail
  } else {
    pred_0 <- rowSums(mod_mat %*% diag(coefs)) + cox_tsrif$frail
  }

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
