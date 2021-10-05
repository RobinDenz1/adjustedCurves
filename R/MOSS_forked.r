##################### MOSS functions with some reworks #########################

# This code was written by Wilson Cai and taken from the official
# github repository: https://github.com/wilsoncai1992/MOSS
# I made a few tiny changes but the main reason this copy is here is that
# CRAN does not allow dependencies on packages that are not on CRAN.
# Since MOSS runs under a GPL-2 license and the adjustedCurves package
# does so as well, this is completely fine legally speaking

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

utils::globalVariables(c("Q1Haz", "G_dC"))

MOSS_hazard <- R6::R6Class("MOSS_hazard",
  public = list(
    A = NULL,
    T_tilde = NULL,
    Delta = NULL,
    density_failure = NULL,
    density_censor = NULL,
    g1W = NULL,
    A_intervene = NULL,

    epsilon = NULL,
    max_num_interation = NULL,
    tmle_tolerance = NULL,
    k_grid = NULL,

    q_best = NULL,
    initialize = function(
      A,
      T_tilde,
      Delta,
      density_failure,
      density_censor,
      g1W,
      A_intervene = NULL,
      k_grid = NULL
    ) {
      self$A <- A
      self$T_tilde <- T_tilde
      self$Delta <- Delta
      self$density_failure <- density_failure
      self$density_censor <- density_censor
      self$g1W <- g1W
      self$A_intervene <- A_intervene

      self$k_grid <- k_grid
      return(self)
    },
    create_dNt = function() {
      dNt <- matrix(0, nrow = length(self$A), ncol = max(self$T_tilde))
      for (i in seq_len(length(self$A))) {
        if (self$Delta[i] == 1) {
          dNt[i, self$T_tilde[i]] <- 1
        }
      }
      return(as.vector(t(dNt)))
    },
    construct_long_data = function(A_intervene, density_failure,
                                   density_censor, which_A="obs") {
      psi_n <- colMeans(density_failure$survival)
      if (which_A == "obs") {
        A <- self$A
      } else {
        A <- self$A_intervene
      }
      eic_fit <- eic$new(
        A = A,
        T_tilde = self$T_tilde,
        Delta = self$Delta,
        density_failure = density_failure,
        density_censor = density_censor,
        g1W = self$g1W,
        psi = psi_n,
        A_intervene = A_intervene
      )
      k_grid <- 1:max(self$T_tilde)
      h_matrix <- list()
      for (k in k_grid) {
        h <- eic_fit$clever_covariate(k = k)
        h_matrix <- c(h_matrix, list(h))
      }
      h_matrix <- do.call(cbind, h_matrix)
      return(h_matrix)
    },
    fit_epsilon = function(method = "l2", clipping = Inf) {
      dNt <- self$create_dNt()
      h_matrix <- self$construct_long_data(
        A_intervene = self$A_intervene,
        density_failure = self$density_failure,
        density_censor = self$density_censor
      )
      if (method == "glm") {
        submodel_fit <- glm.fit(
          x = h_matrix,
          y = dNt,
          family = binomial(),
          offset = logit(as.vector(t(self$density_failure$hazard))),
          intercept = FALSE
        )
        epsilon_n <- submodel_fit$coefficients
        l2_norm <- sqrt(sum(epsilon_n ^ 2))
        if (l2_norm >= clipping) {
          # clipping the step size
          epsilon_n <- epsilon_n / l2_norm * clipping
        }
      }
      if (method %in% c("l2", "l1")) {
        epsilon_n <- tryCatch({
          if (method == "l2") {alpha <- 0
                               norm_func <- norm_l2
                               lambda.min.ratio <- 1e-2}
          if (method == "l1") {alpha <- 1
                               norm_func <- norm_l1
                               lambda.min.ratio <- 9e-1}
          ind <- 1
          while (ind == 1) {
            enet_fit <- glmnet::glmnet(
              x = h_matrix,
              y = dNt,
              offset = logit(as.vector(t(self$density_failure$hazard))),
              family = "binomial",
              alpha = alpha,
              standardize = FALSE,
              intercept = FALSE,
              lambda.min.ratio = lambda.min.ratio,
              # nlambda = 2e2
              nlambda = 1e2
            )
            norms <- apply(enet_fit$beta, 2, norm_func)
            ind <- max(which(norms <= clipping))
            # browser()
            if (ind > 1) break
            # lambda.min.ratio <- (lambda.min.ratio + 1) / 2
            lambda.min.ratio <- sort(enet_fit$lambda, decreasing = TRUE)[2] /
                                     max(enet_fit$lambda)
          }
          epsilon_n <- enet_fit$beta[, ind]
        }, error = function(e) {
          # if error, epsilon = 0
          return(rep(0, ncol(h_matrix)))
        })
      }
      h_matrix_update <- self$construct_long_data(
        A_intervene = self$A_intervene,
        density_failure = self$density_failure,
        density_censor = self$density_censor,
        which_A = self$A_intervene
      )
      hazard_new <- expit(
        logit(as.vector(t(self$density_failure$hazard))) +
          as.vector(h_matrix_update %*% epsilon_n)
      )
      hazard_new <- matrix(
        hazard_new,
        nrow = length(self$A),
        ncol = max(self$T_tilde),
        byrow = TRUE
      )
      # the new hazard for failure
      return(
        survival_curve$new(
          t = 1:max(self$T_tilde), hazard = hazard_new
        )$hazard_to_survival()
      )
    },
    compute_mean_eic = function(psi_n, k_grid) {
      eic_fit <- eic$new(
        A = self$A,
        T_tilde = self$T_tilde,
        Delta = self$Delta,
        density_failure = self$density_failure,
        density_censor = self$density_censor,
        g1W = self$g1W,
        psi = psi_n,
        A_intervene = self$A_intervene
      )$all_t(k_grid = k_grid)
      mean_eic <- colMeans(eic_fit)
      return(mean_eic)
    },
    iterate_onestep = function(
      method = "l2",
      epsilon = 1e0,
      max_num_interation = 1e2,
      tmle_tolerance = NULL,
      verbose = FALSE
    ) {
      self$epsilon <- epsilon
      self$max_num_interation <- max_num_interation
      if (is.null(tmle_tolerance)) {
        self$tmle_tolerance <- 1 / self$density_failure$n()
      } else {
        self$tmle_tolerance <- tmle_tolerance
      }
      k_grid <- 1:max(self$T_tilde)

      psi_n <- colMeans(self$density_failure$survival)
      mean_eic <- self$compute_mean_eic(psi_n = psi_n, k_grid = k_grid)

      num_iteration <- 0

      mean_eic_inner_prod_prev <- abs(sqrt(sum(mean_eic ^ 2)))
      mean_eic_inner_prod_current <- mean_eic_inner_prod_prev
      mean_eic_inner_prod_best <- sqrt(sum(mean_eic ^ 2))
      self$q_best <- self$density_failure$clone(deep = TRUE)
      to_iterate <- TRUE
      if (is.infinite(mean_eic_inner_prod_current) |
          is.na(mean_eic_inner_prod_current)) {
        to_iterate <- FALSE
      }

      while (
        mean_eic_inner_prod_current >= self$tmle_tolerance * sqrt(max(k_grid)) &
        to_iterate
      ) {
        if (verbose) {
          df_debug <- data.frame(num_iteration,
                                 sqrt(sum(mean_eic ^ 2)),
                                 mean(psi_n))
          colnames(df_debug) <- NULL
          print(df_debug)
        }
        # update
        self$density_failure <- self$fit_epsilon(
          method = method, clipping = self$epsilon
        )

        psi_n <- colMeans(self$density_failure$survival)
        mean_eic <- self$compute_mean_eic(psi_n = psi_n, k_grid = k_grid)
        # new stopping
        mean_eic_inner_prod_prev <- mean_eic_inner_prod_current
        mean_eic_inner_prod_current <- abs(sqrt(sum(mean_eic ^ 2)))
        num_iteration <- num_iteration + 1
        if (is.infinite(mean_eic_inner_prod_current) |
            is.na(mean_eic_inner_prod_current)) {
          warning("stopping criteria diverged. Reporting best result so far.")
          break
        }
        if (mean_eic_inner_prod_current < mean_eic_inner_prod_best) {
          # the update caused PnEIC to beat the current best
          # update our best candidate
          self$q_best <- self$density_failure$clone(deep = TRUE)
          mean_eic_inner_prod_best <- mean_eic_inner_prod_current
        }
        if (num_iteration == self$max_num_interation) {
          warning("Max number of iteration reached, stop TMLE")
          break
        }
      }
      # always output the best candidate for final result
      self$density_failure <- self$q_best
      psi_n <- colMeans(self$density_failure$survival)
      if (verbose) {
        message(paste(
          "Pn(EIC)=",
          formatC(mean_eic_inner_prod_best, format = "e", digits = 2),
          "Psi=",
          formatC(mean(psi_n), format = "e", digits = 2)
        ))
      }
      return(psi_n)
    }
  )
)

compute_q <- function(corr, B = 1e3, alpha = 0.05) {
  dim <- nrow(corr)
  z <- apply(
    abs(MASS::mvrnorm(B, mu = rep(0, dim), Sigma = corr)), 1, max
  )
  return(as.numeric(stats::quantile(z, 1 - alpha)))
}

## re-define and rename "compute_simultaneous_ci" to allow different
# alpha levels
compute_se_moss <- function(eic_fit, alpha) {
  # compute the value to +- around the Psi_n
  n <- nrow(eic_fit)
  sigma_squared <- stats::cov(eic_fit)
  sigma <- stats::cor(eic_fit)
  # impute when the variance are zero
  sigma_squared[is.na(sigma_squared)] <- 1e-10
  sigma[is.na(sigma)] <- 1e-10

  variance_marginal <- diag(sigma_squared)
  q <- compute_q(corr=sigma, B=1e3, alpha=alpha)
  return(sqrt(variance_marginal) / sqrt(n) * q)
}


initial_sl_fit <- function(
  T_tilde,
  Delta,
  A,
  W,
  t_max,
  sl_failure = c("SL.glm"),
  sl_censoring = c("SL.glm"),
  sl_treatment = c("SL.glm"),
  gtol = 1e-3
) {
  # convert dictionary of variable names
  ftime <- T_tilde
  ftype <- Delta
  trt <- A
  adjustVars <- W
  t_0 <- t_max
  trtOfInterest <- 0:1
  SL.ftime <- sl_failure
  SL.ctime <- sl_censoring
  SL.trt <- sl_treatment

  adjustVars <- data.frame(adjustVars)
  ftypeOfInterest <- unique(ftype)
  n <- length(ftime)
  id <- seq_len(n)
  dat <- data.frame(id = id, ftime = ftime, ftype = ftype, trt = trt)
  if (!is.null(adjustVars)) dat <- cbind(dat, adjustVars)

  nJ <- length(ftypeOfInterest)
  allJ <- sort(unique(ftype[ftype != 0]))
  ofInterestJ <- sort(ftypeOfInterest)

  # calculate number of groups
  ntrt <- length(trtOfInterest)
  uniqtrt <- sort(trtOfInterest)

  # estimate trt probabilities
  trtOut <- survtmle::estimateTreatment(
    dat = dat,
    ntrt = ntrt,
    uniqtrt = uniqtrt,
    adjustVars = adjustVars,
    SL.trt = SL.trt,
    returnModels = TRUE,
    gtol = gtol
  )
  dat <- trtOut$dat
  trtMod <- trtOut$trtMod

  # make long version of data sets needed for estimation and prediction
  dataList <- survtmle::makeDataList(
    dat = dat, J = allJ, ntrt = ntrt, uniqtrt = uniqtrt, t0 = t_0,
    bounds = NULL
  )
  # estimate censoring
  # when there is almost no censoring, the classification will fail;
  # we manually input the conditional survival for the censoring
  censOut <- tryCatch({
    survtmle::estimateCensoring(
      dataList = dataList,
      ntrt = ntrt,
      uniqtrt = uniqtrt,
      t0 = t_0,
      verbose = FALSE,
      adjustVars = adjustVars,
      SL.ctime = SL.ctime,
      glm.family = "binomial",
      returnModels = TRUE,
      gtol = gtol
    )
  },
  error = function(cond) {
    message("censoring sl error")
    NULL
  }
  )
  if (is.null(censOut)) {
    censOut <- list()
    censOut$dataList <- dataList
    censOut$dataList$obs[, "G_dC"] <- 1
    censOut$dataList$'0'[, "G_dC"] <- 1
    censOut$dataList$'1'[, "G_dC"] <- 1
    is_sl_censoring_converge <- FALSE
    dataList <- censOut$dataList
  } else {
    dataList <- censOut$dataList
    ctimeMod <- censOut$ctimeMod
    is_sl_censoring_converge <- TRUE
  }

  # estimate cause specific hazards
  estOut <- survtmle::estimateHazards(
    dataList = dataList,
    J = allJ,
    verbose = FALSE,
    bounds = NULL,
    adjustVars = adjustVars,
    SL.ftime = SL.ftime,
    glm.family = "binomial",
    returnModels = TRUE
  )
  dataList <- estOut$dataList
  ftimeMod <- estOut$ftimeMod
  # check for convergence
  suppressWarnings(
    if (all(dataList[[1]] == "convergence failure")) {
      return("estimation convergence failure")
    }
  )

  # extract g
  g_1 <- dat$g_1
  g_0 <- dat$g_0

  # extract hazard
  d1 <- dataList$`1`
  d0 <- dataList$`0`

  haz1 <- d1[, c("id", "t", "Q1Haz")]
  haz1 <- tidyr::spread(haz1, t, Q1Haz)
  haz1$id <- NULL # remove the id column

  haz0 <- d0[, c("id", "t", "Q1Haz")]
  haz0 <- tidyr::spread(haz0, t, Q1Haz)
  haz0$id <- NULL # remove the id column

  # extract S_{Ac}
  S_Ac_1 <- d1[, c("id", "t", "G_dC")]
  S_Ac_1 <- tidyr::spread(S_Ac_1, t, G_dC)
  S_Ac_1 <- S_Ac_1[, -1] # remove the id column

  S_Ac_0 <- d0[, c("id", "t", "G_dC")]
  S_Ac_0 <- tidyr::spread(S_Ac_0, t, G_dC)
  S_Ac_0 <- S_Ac_0[, -1] # remove the id column

  density_failure_1 <- survival_curve$new(
    t = seq(range(ftime)[1], range(ftime)[2]), hazard = haz1
  )
  density_failure_0 <- survival_curve$new(
    t = seq(range(ftime)[1], range(ftime)[2]), hazard = haz0
  )
  density_censor_1 <- survival_curve$new(
    t = seq(range(ftime)[1], range(ftime)[2]), survival = S_Ac_1
  )
  density_censor_0 <- survival_curve$new(
    t = seq(range(ftime)[1], range(ftime)[2]), survival = S_Ac_0
  )
  return(list(
    density_failure_1 = density_failure_1,
    density_failure_0 = density_failure_0,
    density_censor_1 = density_censor_1,
    density_censor_0 = density_censor_0,
    g1W = g_1[, 1]
  ))
}

survival_curve <- R6::R6Class("survival_curve",
  public = list(
    t = NULL,
    hazard = NULL,
    # P( T >= t)
    survival = NULL,
    pdf = NULL,
    initialize = function(t, hazard = NULL, survival = NULL, pdf = NULL) {
      # only supports integer grid
      from_hazard <- !is.null(hazard)
      from_survival <- !is.null(survival)
      from_pdf <- !is.null(pdf)
      if (from_hazard + from_survival + from_pdf > 1) {
        stop("cannot construct from both")
      }
      if (!all.equal(t, seq(range(t)[1], range(t)[2]))) {
        stop("t is not integer without gap")
      }
      self$t <- t
      if (from_hazard) {
        # message("construct from hazard")
        if ("data.frame" %in% class(hazard)) hazard <- as.matrix(hazard)
        if ("numeric" %in% class(hazard)) hazard <- matrix(hazard, nrow = 1)
        self$hazard <- hazard
      }
      if (from_survival) {
        # message("construct from survival")
        if ("data.frame" %in% class(survival)) survival <- as.matrix(survival)
        if ("numeric" %in% class(survival)) survival <- matrix(survival,
                                                               nrow = 1)
        self$survival <- survival
      }
      if (from_pdf) {
        # message("construct from pdf")
        if ("data.frame" %in% class(pdf)) pdf <- as.matrix(pdf)
        if ("numeric" %in% class(pdf)) pdf <- matrix(pdf, nrow = 1)
        self$pdf <- pdf
      }
    },
    n = function() {
      n1 <- nrow(self$hazard)
      n2 <- nrow(self$survival)

      return(ifelse(is.null(n1), n2, n1))
    },
    hazard_to_survival = function() {
      # working
      self$survival <- matrix(NA, nrow = self$n(), ncol = max(self$t))
      for (i in 1:self$n()) {
        hazard_here <- c(0, self$hazard[i, ])
        hazard_here <- hazard_here[-length(hazard_here)]
        self$survival[i, ] <- cumprod(1 - hazard_here)
      }
      return(self)
    },
    hazard_to_pdf = function() {
      self$hazard_to_survival()
      # not good using the theory formula
      # self$pdf <- self$hazard * self$survival
      self$survival_to_pdf()
      return(self)
    },
    pdf_to_survival = function() {
      pdf2 <- cbind(0, self$pdf)
      pdf2 <- pdf2[, -ncol(pdf2)]
      # transpose: so that each row is one curve
      self$survival <- 1 - t(apply(pdf2, 1, cumsum))
    },
    pdf_to_hazard = function() {
      self$pdf_to_survival()
      self$hazard <- self$pdf / self$survival
    },
    survival_to_pdf = function() {
      self$pdf <- matrix(NA, nrow = self$n(), ncol = max(self$t))
      for (i in 1:self$n()) {
        self$pdf[i, ] <- c(- diff(self$survival[i, ]), 0)
      }
      return(self)
    },
    survival_to_hazard = function() {
      self$survival_to_pdf()
      self$hazard <- self$pdf / self$survival
      return(self)
    },
    display = function(type, W = NULL) {
      if (is.null(W)) {
        df <- data.frame(t = rep(self$t, self$n()))
      } else {
        if (class(W) != "numeric") stop("W only be univariate vector")
        if (length(W) != self$n()) stop("W length not correct")
        # the first Tmax rows are for the first subject
        df <- data.frame(
          t = rep(self$t, self$n()),
          W = rep(W, each = length(self$t))
        )
      }
      if (type == "survival") {
        df$s <- as.vector(t(self$survival))
        if (!is.null(W)) {
          gg <- ggplot(df, aes(x = t, y = round(W, digits = 1), z = s)) +
            geom_raster(aes(fill = s), interpolate = TRUE) +
            xlim(c(1, max(self$t))) +
            ylab("W") +
            theme_bw()
        } else {
          gg <- ggplot(df, aes(x = t, y = s)) +
            geom_line() +
            theme_bw()
          # ylim(c(-.1, 1.1))
        }
      }
      if (type == "hazard") {
        df$hazard <- as.vector(t(self$hazard))
        if (!is.null(W)) {
          gg <- ggplot(df, aes(x = t, y = round(W, digits = 1), z = hazard)) +
            geom_raster(aes(fill = hazard), interpolate = TRUE) +
            xlim(c(1, max(self$t))) +
            ylab("W") +
            theme_bw()
        } else {
          gg <- ggplot(df, aes(x = t, y = hazard)) +
            geom_line() +
            theme_bw()
        }
      }
      if (type == "pdf") {}
      return(gg)
    },
    create_ggplot_df = function(W = NULL) {
      if (is.null(W)) {
        # only for marginal survival curve
        return(data.frame(t = self$t, s = as.numeric(self$survival)))
      } else {
        if (class(W) != "numeric") stop("W only be univariate vector")
        if (length(W) != self$n()) stop("W length not correct")
        # the first Tmax rows are for the first subject
        df <- data.frame(
          t = rep(self$t, self$n()),
          W = rep(W, each = length(self$t)),
          s = as.vector(t(self$survival))
        )
        return(df)
      }
    },
    ci = function(
      A,
      T_tilde,
      Delta,
      density_failure,
      density_censor,
      g1W,
      psi_n,
      A_intervene,
      alpha = 0.05
    ) {
      eic_fit <- eic$new(
        A = A,
        T_tilde = T_tilde,
        Delta = Delta,
        density_failure = density_failure,
        density_censor = density_censor,
        g1W = g1W,
        psi = psi_n,
        A_intervene = A_intervene
      )$all_t(k_grid = self$t)
      sigma <- apply(eic_fit, 2, sd)
      lower <- psi_n - sigma * 1.96
      upper <- psi_n + sigma * 1.96
      return(data.frame(t = self$t, lower = lower, upper = upper))
    }
  )
)

eic <- R6::R6Class("eic",
  public = list(
    A = NULL,
    T_tilde = NULL,
    Delta = NULL,
    density_failure = NULL,
    density_censor = NULL,
    g1W = NULL,
    psi = NULL,
    A_intervene = NULL,
    initialize = function(
      A, T_tilde, Delta, density_failure, density_censor, g1W, psi, A_intervene
    ) {
      self$A <- A
      self$T_tilde <- T_tilde
      self$Delta <- Delta
      self$density_failure <- density_failure
      self$density_censor <- density_censor
      self$g1W <- g1W
      self$psi <- psi
      self$A_intervene <- A_intervene
      return(self)
    },
    one_t = function(k) {
      if (self$A_intervene == 1) g <- self$g1W  else g <- 1 - self$g1W
      part1_sum <- rep(0, length(g))
      for (t in 1:k) {
        h <- -as.numeric(self$A == self$A_intervene) / g /
          self$density_censor$survival[, t] *
          self$density_failure$survival[, k] /
          self$density_failure$survival[, t]
        part1 <-  h * (
          as.numeric(self$T_tilde == t & self$Delta == 1) -
            as.numeric(self$T_tilde >= t) * self$density_failure$hazard[, t]
        )
        part1_sum <- part1_sum + part1
      }
      part2 <- self$density_failure$survival[, k] - self$psi[k]
      return(part1_sum + part2)
    },
    all_t = function(k_grid) {
      # naive way to compute for all t
      eic_all <- list()
      for (k in k_grid) {
        eic_all <- c(eic_all, list(self$one_t(k = k)))
      }
      eic_all <- do.call(cbind, eic_all)
      return(eic_all)
    },
    clever_covariate = function(k) {
      if (self$A_intervene == 1) g <- self$g1W  else g <- 1 - self$g1W
      h_list <- list()
      for (t in 1:max(self$T_tilde)) {
        if (t > k) {
          # clever covariate is zero beyond
          h <- rep(0, length(g))
        } else {
          h <- -as.numeric(self$A == self$A_intervene) / g /
            self$density_censor$survival[, t] *
            self$density_failure$survival[, k] /
            self$density_failure$survival[, t]
        }
        h_list <- c(h_list, list(h))
      }
      # the first row is 1 ~ t_max for the first subject
      h_list <- do.call(cbind, h_list)
      # the first 1 ~ t_max element is for the first subject
      return(as.vector(t(h_list)))
    }
  )
)

expit <- function(x) exp(x) / (1 + exp(x))

logit <- function(x) log(x) - log(1 - x)

norm_l2 <- function(beta) sqrt(sum(beta ^ 2))

norm_l1 <- function(beta) sum(abs(beta))
