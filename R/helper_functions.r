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

## estimate iptw weights
get_iptw_weights <- function(data, treatment_model, weight_method,
                             variable, stabilize=T, ...) {

  # using WeightIt
  if (inherits(treatment_model, "formula")) {
    args <- list(formula=treatment_model, data=data,
                 method=weight_method,
                 estimand="ATE")
    weights <- do.call(WeightIt::weightit, c(args, ...))$weights
  # using a logistic regression model
  } else if (inherits(treatment_model, "glm")) {
    ps <- stats::predict.glm(treatment_model, newdata=data, type="response")
    weights <- ifelse(data[, variable]==1, 1/ps, 1/(1-ps))
  # using a multinomial logistic regression model
  } else if (inherits(treatment_model, "multinom")) {
    predict.multinom <- utils::getFromNamespace("predict.multinom", "nnet")
    ps <- predict.multinom(treatment_model, newdata=data, type="probs")

    weights <- rep(0, nrow(data))
    for (i in levels(data[,variable])) {
      weights[data[,variable] == i] <- 1/ps[data[,variable] == i, i]
    }
  # nothing else allowed
  } else {
    stop("Unsuported input: '", class(treatment_model), "'. See documentation.")
  }

  if (stabilize) {
    weights <- weights * length(weights) / sum(weights)
  }

  return(weights)
}

## a stat needed to calculate the variance of the iptw_km method
calc_Mj <- function(t, data, ev_time, ps_score) {

  relevant_ps <- ps_score[data[,ev_time] >= t]
  Mj <- ((sum(1/relevant_ps))^2) / (sum((1/relevant_ps)^2))
  return(Mj)

}

## calculate the variance of iptw_km estimates
calc_iptw_km_var <- function(t, adj_km) {

  rel_t <- adj_km[adj_km$time <= t,]
  rel_t$in_sum <- (1-rel_t$s_j) / (rel_t$Mj * rel_t$s_j)
  var_iptw_km <- rel_t$surv[rel_t$time==t]^2 * sum(rel_t$in_sum)
  return(var_iptw_km)

}

## Computes the standard error of a weighted mean using one of
## four possible approximations
weighted.var.se <- function(x, w, se_method, na.rm=F) {

  if (na.rm) {
    miss_ind <- !is.na(x)
    w <- w[miss_ind]
    x <- x[miss_ind]
  }

  n <- length(x)
  mean_Xw <- stats::weighted.mean(x=x, w=w, na.rm=na.rm)

  ## Miller (1977)
  if (se_method=="miller") {
    se <- 1/n * (1/sum(w)) * sum(w * (x - mean_Xw)^2)
  ## Galloway et al. (1984)
  } else if (se_method=="galloway") {
    se <- (n/(sum(w)^2)) * ( (n*sum(w^2 * x^2) - sum(w*x)^2) / (n*(n-1)) )
  ## Cochrane (1977)
  } else if (se_method=="cochrane") {
    mean_W <- mean(w)
    se <- (n/((n-1)*sum(w)^2))*(sum((w*x - mean_W*mean_Xw)^2)
                                - 2*mean_Xw*sum((w-mean_W)*(w*x-mean_W*mean_Xw))
                                + mean_Xw^2*sum((w-mean_W)^2))
  ## just use a classic weighted version (Hmisc)
  } else if (se_method=="simple") {
    se <- (sum(w * (x - mean_Xw)^2) / (sum(w) - 1)) / n
  }
  return(se)
}

## function to get predicted values from any kind of geese object,
## for multiple time points
geese_predictions <- function(geese_mod, Sdata, times, n) {

  # full model matrix and betas
  mod_mat <- stats::model.matrix(geese_mod$formula,
                                 data=stats::model.frame(geese_mod$formula, Sdata,
                                             na.action=stats::na.pass))
  betas <- geese_mod$beta

  apply_betas <- function(x, betas) {
    return(sum(x * betas))
  }

  pred_mat <- matrix(nrow=n, ncol=length(times))
  for (i in 1:length(times)) {
    # take only relevant portion (at time t) of model matrix
    mod_mat_t <- mod_mat[Sdata$vtime==times[i],]
    # apply coefficients
    preds <- apply(X=mod_mat_t, MARGIN=1, FUN=apply_betas, betas=betas)
    pred_mat[,i] <- preds
  }
  colnames(pred_mat) <- paste0("t.", times)
  return(pred_mat)
}

## calculate CI from sd
confint_surv <- function(surv, se, conf_level, conf_type="plain") {
  # critical value
  z_val <- stats::qnorm(1-((1-conf_level)/2))

  if (conf_type=="plain") {
    error <- z_val * se
    left <- surv - error
    right <- surv + error
  } else if (conf_type=="log") {
    error <- z_val * se
    left <- surv * exp(-error)
    right <- surv * exp(error)
  } else if (conf_type=="log-log") {
    xx <- ifelse(surv==0 | surv==1, NA, surv)
    se2 <- z_val * se/log(xx)
    left <- exp(-exp(log(-log(xx)) - se2))
    right <- exp(-exp(log(-log(xx)) + se2))
  }

  return(list(left=left, right=right))
}

## takes a value x at which to read from the step function
## and step function data from which to read it
read_from_step_function <- function(x, step_data, est="surv") {

  # no extrapolation
  if (x > max(step_data$time)) {
    return(NA)
  }

  # otherwise get value
  check <- step_data[which(step_data$time <= x),]
  if (nrow(check)==0) {
    if (est=="surv") {
      val <- 1
    } else {
      val <- 0
    }
  } else {
    val <- check[,est][which(check$time==max(check$time))][1]
  }
  return(val)
}

## calculate difference between two step functions
## according to some transformation function
exact_stepfun_difference <- function(adjsurv, times, est="surv") {

  levs <- unique(adjsurv$group)
  adjsurv_0 <- adjsurv[which(adjsurv$group==levs[1]),]
  adjsurv_1 <- adjsurv[which(adjsurv$group==levs[2]),]

  if (nrow(adjsurv_0) == nrow(adjsurv_1)) {

    if (all(adjsurv_0$time == adjsurv_1$time)) {
      surv_0 <- adjsurv_0$surv
      surv_1 <- adjsurv_1$surv
    } else {
      surv_0 <- sapply(times, read_from_step_function, step_data=adjsurv_0,
                       est=est)
      surv_1 <- sapply(times, read_from_step_function, step_data=adjsurv_1,
                       est=est)
    }

  } else {
    surv_0 <- sapply(times, read_from_step_function, step_data=adjsurv_0,
                     est=est)
    surv_1 <- sapply(times, read_from_step_function, step_data=adjsurv_1,
                     est=est)
  }

  surv_diff <- surv_1 - surv_0

  diff_dat <- data.frame(time=times)
  diff_dat[,est] <- surv_diff

  return(diff_dat)
}

## calculate exact integral under step function
# 'stepfun' needs to be a data.frame with columns 'time' and 'surv',
# sorted by time with no duplicates in time
exact_stepfun_integral <- function(stepfun, from, to, est="surv") {

  # constrain step function end
  latest <- read_from_step_function(to, step_data=stepfun, est=est)
  stepfun <- stepfun[stepfun$time <= to,]

  if (!to %in% stepfun$time) {
    temp <- data.frame(time=to)
    temp[,est] <- latest
    stepfun <- rbind(stepfun, temp)
  }

  # constrain step function beginning
  if (from != 0) {
    earliest <- read_from_step_function(from, step_data=stepfun, est=est)
    stepfun <- stepfun[stepfun$time >= from,]

    if (!from %in% stepfun$time) {
      temp <- data.frame(time=from)
      temp[,est] <- earliest
      stepfun <- rbind(temp, stepfun)
    }

  }

  # when there are unknown survival times, return NA
  # only neccessary when the last time is NA, since my
  # algorithm technically still works in that case, but makes no sense
  if (anyNA(stepfun)) {
    return(NA)
  }

  # calculate exact integral
  integral <- 0
  for (i in 1:(length(stepfun$time)-1)) {
    x1 <- stepfun$time[i]
    x2 <- stepfun$time[i+1]
    y <- stepfun[,est][i]
    rect_area <- (x2 - x1) * y
    integral <- integral + rect_area
  }
  return(integral)
}

## function to change plotdata in iptw and standard methods when
## custom points in time are supplied
specific_times <- function(plotdata, times, cif=F) {

  levs <- unique(plotdata$group)
  new_plotdata <- vector(mode="list", length=length(levs))
  for (i in 1:length(levs)) {

    if (cif) {
      new_est <- sapply(times, read_from_step_function, est="cif",
                        step_data=plotdata[which(plotdata$group==levs[i]),])
      new_dat <- data.frame(time=times, cif=new_est, group=levs[i])
    } else {

      new_est <- sapply(times, read_from_step_function, est="surv",
                        step_data=plotdata[which(plotdata$group==levs[i]),])
      new_dat <- data.frame(time=times, surv=new_est, group=levs[i])
    }

    if ("se" %in% colnames(plotdata)) {
      # read from curve using custom function
      new_se <- sapply(times, read_from_step_function, est="se",
                       step_data=plotdata[which(plotdata$group==levs[i]),])
      new_ci_lower <- sapply(times, read_from_step_function, est="ci_lower",
                             step_data=plotdata[which(plotdata$group==levs[i]),])
      new_ci_upper <- sapply(times, read_from_step_function, est="ci_upper",
                             step_data=plotdata[which(plotdata$group==levs[i]),])
      # add to output in same order
      new_dat$se <- new_se
      new_dat$ci_lower <- new_ci_lower
      new_dat$ci_upper <- new_ci_upper

    }
    new_plotdata[[i]] <- new_dat
  }
  new_plotdata <- as.data.frame(dplyr::bind_rows(new_plotdata))
  return(new_plotdata)
}

## used to combine output from foreach
multi_result_class <- function(boot_data=NULL, boot_data_same_t=NULL) {
  me <- list(
    boot_data = boot_data,
    boot_data_same_t = boot_data_same_t
  )

  class(me) <- append(class(me),"multiResultClass")
  return(me)
}

## redefine 'timepoints' from survtmle to fix a bug in there
survtmle.timepoints <- function(object, times, returnModels=FALSE,
                                SL.trt, SL.ctime, SL.ftime,
                                glm.trt, glm.ctime, glm.ftime) {

  if (is.null(object$trtMod)) {
    stop("object must have returnModels = TRUE")
  }

  callList <- as.list(object$call)[-1]
  cglm <- any(class(object$ctimeMod) %in% c("glm", "speedglm")) |
    any(class(object$ctimeMod) == "noCens")

  tglm <- any(class(object$trtMod) %in% c("glm", "speedglm"))
  ftglm <- ifelse(callList$method == "hazard",
                  any(class(object$ftimeMod[[1]]) %in% c(
                    "glm",
                    "speedglm"
                  )), FALSE
  )

  myOpts <- c(
    "t0", "returnModels",
    ifelse(cglm, "glm.ctime", "SL.ctime"),
    ifelse(tglm, "glm.trt", "SL.trt")
  )
  if (callList$method == "hazard") {
    myOpts <- c(myOpts, ifelse(ftglm, "glm.ftime", "SL.ftime"))
  }
  funOpts <- callList[-which(names(callList) %in% myOpts)]

  funOpts$returnModels <- returnModels
  # used glm for censoring?
  if (cglm) {
    funOpts$glm.ctime <- object$ctimeMod
    funOpts$SL.ctime <- NULL
  } else {
    funOpts$SL.ctime <- object$ctimeMod
  }
  # used glm for trt?
  if (tglm) {
    funOpts$glm.trt <- object$trtMod
  } else {
    funOpts$SL.trt <- object$trtMod
  }
  # used glm for ftime
  if (ftglm & callList$method == "hazard") {
    funOpts$glm.ftime <- object$ftimeMod
  } else if (!ftglm & callList$method == "hazard") {
    funOpts$SL.ftime <- object$ftimeMod
  }
  # add in failure times, types, trt, and adjust
  funOpts$ftime <- object$ftime
  funOpts$ftype <- object$ftype
  funOpts$trt <- object$trt
  funOpts$adjustVars <- object$adjustVars

  outList <- vector(mode = "list", length = length(times))
  ct <- 0
  for (i in times) {
    ct <- ct + 1
    funOpts$t0 <- i
    if (all(object$ftime[object$ftype > 0] > i)) {
      outList[[ct]] <- list(
        est = rep(0, length(object$est)),
        var = matrix(
          NA,
          nrow = length(object$est),
          ncol = length(object$est)
        )
      )
    } else {
      if (i != object$t0) {
        outList[[ct]] <- do.call("survtmle", args = funOpts)
      } else {
        outList[[ct]] <- object
      }
    }
  }
  names(outList) <- paste0("t", times)
  class(outList) <- "tp.survtmle"
  return(outList)
}

## calculate pseudo values for the survival function
## using either the standard way or with dependent censoring
calc_pseudo_surv <- function(data, ev_time, event, times, censoring_vars,
                             ipcw.method) {

  # standard pseudo-values, no dependent censoring
  if (is.null(censoring_vars)) {

    # estimate pseudo observations
    hist_formula <- stats::as.formula(paste("prodlim::Hist(", ev_time, ", ",
                                            event, ") ~ 1"))
    pseudo <- prodlim::jackknife(prodlim::prodlim(hist_formula, data=data),
                                 times=times)

  } else {

    requireNamespace("eventglm")
    pseudo_aareg <- utils::getFromNamespace("pseudo_aareg", "eventglm")

    cens_formula <- stats::as.formula(paste0("~ ", paste(censoring_vars,
                                                         collapse=" + ")))
    pseudo_formula <- stats::as.formula(paste0("survival::Surv(", ev_time,
                                               ", ", event, ") ~ 1"))

    pseudo <- sapply(times, FUN=pseudo_aareg,
                     formula=pseudo_formula,
                     cause=1,
                     data=data,
                     type="survival",
                     formula.censoring=cens_formula,
                     ipcw.method=ipcw.method)

  }
  return(pseudo)
}

## function to perform direct standardization using a linear model
## given a vector of pseudo values
lm_direct <- function(x, glm_formula, data, levs, variable,
                      outcome_vars, conf_level) {

  # fit linear regression
  data$yi <- x
  mod <- stats::glm(formula=glm_formula, family="gaussian", data=data)

  # do direct standardization
  # since it's a simple linear model, we can actually use
  # the 'average of covariate' method, which also allows us to
  # get confidence intervals for the estimate directly from the
  # predict.lm function (uses the delta method internally)
  newdata <- as.data.frame(t(apply(data[,outcome_vars], 2, mean)))
  out <- list()
  for (i in 1:length(levs)) {
    newdata[,variable] <- levs[i]
    pred <- stats::predict.lm(mod, newdata=newdata, se.fit=T,
                              interval="confidence", level=conf_level)
    pred_dat <- as.data.frame(pred$fit)
    pred_dat$group <- levs[i]
    pred_dat$se <- pred$se.fit
    out[[i]] <- pred_dat
  }
  out <- dplyr::bind_rows(out)
  return(out)
}

## function to model treatment assignment using SuperLearner
get_SL_ps_score <- function(data, variable, treatment_vars, SL.trt,
                            cv_folds) {

  # build prediction model
  form <- paste0("~", paste0(treatment_vars, collapse=" + "))
  SL_dat <- as.data.frame(stats::model.matrix(stats::as.formula(form),
                                              data=data[,treatment_vars]))
  SL_dat$`(Intercept)` <- NULL
  sl_fit <- SuperLearner::SuperLearner(Y=data[,variable], X=SL_dat, newX=NULL,
                                       family="binomial", SL.library=SL.trt,
                                       control=list(saveFitLibrary=TRUE),
                                       cvControl=list(V=cv_folds))

  # estimate propensity score
  ps_score <- SuperLearner::predict.SuperLearner(sl_fit,
                                                 newdata=data,
                                                 X=SL_dat,
                                                 Y=data[,variable])$pred

  return(ps_score)
}

## calculating the variance of tmle_pseudo point estimates
tmle_pseudo_var <- function(Y, A, ps_score, Qstar, n) {

  QA <- ifelse(A==0, Qstar[,1], Qstar[,2])

  mu_0 <- mean(Qstar[,1])
  mu_1 <- mean(Qstar[,2])

  g0W <- 1 - ps_score
  g1W <- ps_score

  # efficient influence curve estimates
  IC_0 <- (1-A)/g0W * (Y - QA) + Qstar[,1] - mu_0
  IC_1 <- A/g1W * (Y - QA) + Qstar[,2] - mu_1

  # variance approximations using the ICs
  return(c(var_0=stats::var(IC_0)/n,
           var_1=stats::var(IC_1)/n))

}
