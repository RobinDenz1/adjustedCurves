####################### Helper functions for R-Package #########################

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
## three possible approximation
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

  current.na.action <- options('na.action')
  options(na.action="na.pass")
  # full model matrix and betas
  mod_mat <- stats::model.matrix(geese_mod$formula, data=Sdata)
  options(na.action=current.na.action[[1]])

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
  if (conf_type=="plain") {
    error <- stats::qnorm(1-((1-conf_level)/2)) * se
    left <- surv - error
    right <- surv + error
  } else if (conf_type=="log") {
    error <- stats::qnorm(1-((1-conf_level)/2)) * se
    left <- surv * exp(-error)
    right <- surv * exp(error)
  }

  return(list(left=left, right=right))
}

## re-define and rename "compute_simultaneous_ci" to allow different alpha levels
compute_se_moss <- function(eic_fit, alpha) {
  # compute the value to +- around the Psi_n
  n <- nrow(eic_fit)
  sigma_squared <- stats::cov(eic_fit)
  sigma <- stats::cor(eic_fit)
  # impute when the variance are zero
  sigma_squared[is.na(sigma_squared)] <- 1e-10
  sigma[is.na(sigma)] <- 1e-10

  variance_marginal <- diag(sigma_squared)
  q <- MOSS::compute_q(corr=sigma, B=1e3, alpha=alpha)
  return(sqrt(variance_marginal) / sqrt(n) * q)
}

## simulate survival time according to Bender et al. (2005)
sim_surv_time <- function(row, betas, dist, lambda, gamma) {
  U <- stats::runif(1, min=0, max=1)
  eff <- sum(row * betas)

  if (dist=="weibull") {
    surv_time <- (-(log(U)/(lambda*exp(eff))))^(1/gamma)
  } else if (dist=="exponential") {
    surv_time <- -(log(U)/(lambda*exp(eff)))
  }
  return(surv_time)
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
    val <- 1
  } else {
    val <- check[,est][which(check$time==max(check$time))][1]
  }
  return(val)
}

## calculate difference between two step functions
## according to some transformation function
exact_stepfun_difference <- function(adjsurv, times) {

  levs <- unique(adjsurv$group)
  adjsurv_0 <- adjsurv[which(adjsurv$group==levs[1]),]
  adjsurv_1 <- adjsurv[which(adjsurv$group==levs[2]),]

  if (nrow(adjsurv_0) == nrow(adjsurv_1)) {

    if (all(adjsurv_0$time == adjsurv_1$time)) {
      surv_0 <- adjsurv_0$surv
      surv_1 <- adjsurv_1$surv
    } else {
      surv_0 <- sapply(times, read_from_step_function, step_data=adjsurv_0)
      surv_1 <- sapply(times, read_from_step_function, step_data=adjsurv_1)
    }

  } else {
    surv_0 <- sapply(times, read_from_step_function, step_data=adjsurv_0)
    surv_1 <- sapply(times, read_from_step_function, step_data=adjsurv_1)
  }

  surv_diff <- surv_1 - surv_0

  diff_dat <- data.frame(time=times, surv=surv_diff)

  return(diff_dat)
}

## calculate exact integral under step function
# 'stepfun' needs to be a data.frame with columns 'time' and 'surv',
# sorted by time with no duplicates in time
exact_stepfun_integral <- function(stepfun, from, to) {

  # constrain step function end
  latest <- read_from_step_function(to, step_data=stepfun)
  stepfun <- stepfun[stepfun$time <= to,]

  if (!to %in% stepfun$time) {
    stepfun <- rbind(stepfun, data.frame(time=to, surv=latest))
  }

  # constrain step function beginning
  if (from != 0) {
    earliest <- read_from_step_function(from, step_data=stepfun)
    stepfun <- stepfun[stepfun$time >= from,]

    if (!from %in% stepfun$time) {
      stepfun <- rbind(data.frame(time=from, surv=earliest), stepfun)
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
    y <- stepfun$surv[i]
    rect_area <- (x2 - x1) * y
    integral <- integral + rect_area
  }
  return(integral)
}

## used to combine output from foreach
multi_result_class <- function(boot_data=NULL, boot_data_same_t=NULL) {
  me <- list(
    boot_data = boot_data,
    boot_data_same_t = boot_data_same_t
  )

  ## Set the name for the class
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}

