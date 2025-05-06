
# Assign to global, to get rid off devtools::check() note
utils::globalVariables(c("gaussian", "id"))

## generalized function that does both
method_aiptw_pseudo <- function(data, variable, ev_time, event, cause,
                                conf_int, conf_level, times, outcome_vars,
                                treatment_model, type_time, spline_df,
                                censoring_vars, ipcw_method, mode) {
  # some constants
  len <- length(times)
  n <- nrow(data)
  group <- data[, variable]

  if (is.numeric(treatment_model)) {
    ps_score <- treatment_model
  } else if (inherits(treatment_model, "glm")) {
    ps_score <- treatment_model$fitted.values
  } else if (inherits(treatment_model, "multinom")) {
    predict.multinom <- utils::getFromNamespace("predict.multinom", "nnet")
    ps_score <- predict.multinom(treatment_model, newdata=data, type="probs")
  }

  # estimate pseudo observations
  pseudo <- calc_pseudo_surv(data=data,
                             ev_time=ev_time,
                             event=event,
                             times=times,
                             censoring_vars=censoring_vars,
                             ipcw.method=ipcw_method,
                             cause=cause)
  # create data for geese
  Sdata <- data.frame(yi=c(pseudo),
                      group=rep(group, len),
                      vtime=rep(times, rep(n, len)),
                      id=rep(1:n, len))

  outcome_vars <- outcome_vars[outcome_vars != variable]
  for (col in outcome_vars) {
    Sdata[, col] <- rep(data[, col], len)
  }

  if (type_time=="factor") {
    Sdata$vtime <- as.factor(Sdata$vtime)
    geese_formula <- paste("yi ~ vtime + ", paste(outcome_vars, collapse=" + "),
                           " + group")
  } else if (type_time=="bs") {
    geese_formula <- paste("yi ~ splines::bs(vtime, df=", spline_df, ") + ",
                           paste(outcome_vars, collapse=" + "), " + group")
  } else if (type_time=="ns") {
    geese_formula <- paste("yi ~ splines::ns(vtime, df=", spline_df, ") + ",
                           paste(outcome_vars, collapse=" + "), " + group")
  }

  # remove rows where pseudo-values are NA for geese
  Sdata_fit <- Sdata[!is.na(Sdata$yi), ]

  # call geese
  geese_mod <- geepack::geese(stats::as.formula(geese_formula), scale.fix=TRUE,
                              data=Sdata_fit, family=gaussian, id=id,
                              jack=FALSE, mean.link="cloglog",
                              corstr="independence")

  # get direct adjustment estimates
  levs <- levels(data[, variable])
  plotdata <- vector(mode="list", length=length(levs))

  for (i in seq_len(length(levs))) {

    Sdata$group <- factor(levs[i], levels=levs)
    pred <- geese_predictions(geese_mod, Sdata, times=times, n=n)

    m <- 1 - exp(-exp(pred))

    # augment estimates using propensity score
    if (length(levs) > 2) {
      group_ind <- ifelse(group==levs[i], 1, 0)
      ps_score_lev <- ps_score[, levs[i]]

      dr <- (pseudo*group_ind-(group_ind-ps_score_lev)*m)/ps_score_lev
      # if binary, use equation from the paper directly
    } else if (i == 1) {
      group <- ifelse(data[, variable]==levs[1], 0, 1)
      dr <- (pseudo*(1-group)+(group-ps_score)*m)/(1-ps_score)
    } else if (i == 2) {
      group <- ifelse(data[, variable]==levs[1], 0, 1)
      dr <- (pseudo*group-(group-ps_score)*m)/ps_score
    }

    surv <- apply(dr, 2, mean, na.rm=TRUE)

    if (conf_int) {

      pseudo_dr_se <- function(x, n, na.rm) {
        sqrt(stats::var(x, na.rm=na.rm) / n)
      }
      surv_se <- apply(dr, 2, pseudo_dr_se, n=n, na.rm=TRUE)

      surv_ci <- confint_surv(surv=surv, se=surv_se, conf_level=conf_level,
                              conf_type="plain")

      plotdata[[i]] <- data.frame(time=times, surv=surv, group=levs[i],
                                  se=surv_se, ci_lower=surv_ci$left,
                                  ci_upper=surv_ci$right)

    } else {
      plotdata[[i]] <- data.frame(time=times, surv=surv, group=levs[i])
    }

  }
  plotdata <- dplyr::bind_rows(plotdata)
  rownames(plotdata) <- NULL

  if (mode=="cif") {
    colnames(plotdata)[colnames(plotdata)=="surv"] <- "cif"
  }

  output <- list(plotdata=plotdata,
                 pseudo_values=pseudo,
                 geese_model=geese_mod)

  if (mode=="cif") {
    class(output) <- "adjustedcif.method"
  } else {
    class(output) <- "adjustedsurv.method"
  }

  return(output)

}

## AIPTW with Pseudo-Values for Survival Curves
#' @export
surv_aiptw_pseudo <- function(data, variable, ev_time, event, conf_int,
                              conf_level=0.95, times, outcome_vars,
                              treatment_model, type_time="factor",
                              spline_df=5, censoring_vars=NULL,
                              ipcw_method="binder") {

  out <- method_aiptw_pseudo(data=data, variable=variable, ev_time=ev_time,
                             event=event, conf_int=conf_int,
                             conf_level=conf_level, times=times,
                             outcome_vars=outcome_vars,
                             treatment_model=treatment_model,
                             type_time=type_time, spline_df=spline_df,
                             censoring_vars=censoring_vars,
                             ipcw_method=ipcw_method, mode="surv", cause=1)
  return(out)
}

## AIPTW with Pseudo-Values for CIF
#' @export
cif_aiptw_pseudo <- function(data, variable, ev_time, event, cause,
                             conf_int, conf_level=0.95, times,
                             outcome_vars, treatment_model,
                             type_time="factor", spline_df=5) {

  out <- method_aiptw_pseudo(data=data, variable=variable, ev_time=ev_time,
                             event=event, cause=cause, conf_int=conf_int,
                             conf_level=conf_level, times=times,
                             outcome_vars=outcome_vars,
                             treatment_model=treatment_model,
                             type_time=type_time, spline_df=spline_df,
                             censoring_vars=NULL, ipcw_method="binder",
                             mode="cif")
  return(out)
}
