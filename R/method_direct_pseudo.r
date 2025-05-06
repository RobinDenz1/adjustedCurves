
# Assign to global, to get rid off devtools::check() note
utils::globalVariables(c("gaussian", "id"))

## generalized function that does both
method_direct_pseudo <- function(data, variable, ev_time, event, cause,
                                 conf_int, conf_level, times, outcome_vars,
                                 type_time, spline_df, censoring_vars,
                                 ipcw_method, mode) {

  # estimate pseudo observations
  pseudo <- calc_pseudo_surv(data=data,
                             ev_time=ev_time,
                             event=event,
                             times=times,
                             censoring_vars=censoring_vars,
                             ipcw.method=ipcw_method,
                             cause=cause)

  # remove "variable" from outcome_vars because it is always included
  outcome_vars <- outcome_vars[outcome_vars!=variable]

  # some constants
  len <- length(times)
  n <- nrow(data)
  group <- data[, variable]

  # create data for geese
  Sdata <- data.frame(yi=1 - c(pseudo),
                      group=rep(group, len),
                      vtime=rep(times, rep(n, len)),
                      id=rep(1:n, len))
  for (col in outcome_vars) {
    Sdata[, col] <- rep(data[, col], len)
  }

  if (type_time=="factor") {
    Sdata$vtime <- as.factor(Sdata$vtime)
    geese_formula <- paste("yi ~ vtime + ",
                           paste(outcome_vars, collapse=" + "),
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

  # initialize outcome df list
  levs <- levels(data[, variable])
  plotdata <- vector(mode="list", length=length(levs))

  # do direct adjustment
  for (i in seq_len(length(levs))) {

    Sdata$group <- factor(levs[i], levels=levs)
    pred <- geese_predictions(geese_mod, Sdata, times=times, n=n)

    m <- exp(-exp(pred))
    surv <- apply(m, 2, mean, na.rm=TRUE)

    plotdata[[i]] <- data.frame(time=times, surv=surv, group=levs[i])

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

## Using Pseudo Observations and Direct Adjustment
#' @export
surv_direct_pseudo <- function(data, variable, ev_time, event,
                               conf_int=FALSE, conf_level=0.95, times,
                               outcome_vars, type_time="factor",
                               spline_df=5, censoring_vars=NULL,
                               ipcw_method="binder") {

  out <- method_direct_pseudo(data=data, variable=variable, ev_time=ev_time,
                              event=event, conf_int=conf_int,
                              conf_level=conf_level, times=times,
                              outcome_vars=outcome_vars, type_time=type_time,
                              spline_df=spline_df, cause=1, mode="surv",
                              censoring_vars=censoring_vars,
                              ipcw_method=ipcw_method)
  return(out)
}

## Using Pseudo Observations and Direct Adjustment
#' @export
cif_direct_pseudo <- function(data, variable, ev_time, event, cause,
                              conf_int=FALSE, conf_level=0.95, times,
                              outcome_vars, type_time="factor", spline_df=5) {

  out <- method_direct_pseudo(data=data, variable=variable, ev_time=ev_time,
                              event=event, conf_int=conf_int,
                              conf_level=conf_level, times=times,
                              outcome_vars=outcome_vars, type_time=type_time,
                              spline_df=spline_df, cause=cause, mode="cif",
                              censoring_vars=NULL, ipcw_method="binder")
  return(out)
}
