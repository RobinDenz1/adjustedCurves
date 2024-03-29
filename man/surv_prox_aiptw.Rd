\name{surv_prox_aiptw}
\alias{surv_prox_aiptw}

\title{
Proximal Augmented Inverse Probability of Treatment Weighted Survival Curve Estimates
}
\description{

This page explains the details of estimating augmented inverse probability of treatment weighted survival curves using a proximal causal inference based method for single event time-to-event data (\code{method="prox_aiptw"} as described by Ying et al. (2022) in the \code{\link[=adjustedsurv]{adjustedsurv}} function). All regular arguments of the \code{adjustedsurv} function can be used. Additionally, the \code{treatment_proxy}, \code{outcome_proxy} and \code{adjust_vars} arguments have to be specified in the \code{adjustedsurv} call. Further arguments specific to this method are listed below.

}

\arguments{
  \item{adjust_vars}{
[\strong{required}] A character vector specifying names of variables in \code{data}. These variables may consist of observed confounders. At least one variable has to be named. Can be numeric, character or factor variables. Corresponds to \eqn{X} (type 1 proxy) in the article by Ying et al. (2022).
}
  \item{treatment_proxy}{
[\strong{required}] A single character string specifying the (numeric) variable that should be used as a treatment proxy. Corresponds to \eqn{Z} (type 3 proxy) in the article by Ying et al. (2022).
}
  \item{outcome_proxy}{
[\strong{required}] A single character string specifying the (numeric) variable that should be used as a outcome proxy. Corresponds to \eqn{W} (type 2 proxy) in the article by Ying et al. (2022).
  }
  \item{optim_method}{
A single character string passed to the \code{method} argument of the \code{optim} function, used internally when fitting the q-confounding bridge function and the h-confounding bridge function. Defaults to \code{"BFGS"}. To pass additional argument to the internal \code{optim} call, see argument \code{optim_control}.
  }
  \item{optim_control}{
A list of named arguments passed to the \code{control} argument of the \code{optim} function, used internally when fitting the q-confounding bridge function and the h-confounding bridge function. Set to \code{list()} to not pass any additional argument (default).
  }
  \item{return_fit}{
Whether to add intermediate results, such as the q-confounding bridge function to the output object. Defaults to \code{TRUE}.
  }
}
\details{

\itemize{
\item{\strong{Type of Adjustment:} Uses the proximal causal inference framework to adjust for measured and unmeasured confounding through the use of both the q-confounding bridge function and the h-confounding bridge function, which is essentially augmented inverse probability of treatment weighting, but using proxies.}
\item{\strong{Doubly-Robust:} Estimates are Doubly-Robust in the sense that only one of the bridge functions has to be correctly specified to achieve unbiased estimates, given that the other relevant assumptions hold.}
\item{\strong{Categorical groups:}\code{variable} may only contain two groups. Must be a factor variable.}
\item{\strong{Approximate Variance:} Calculations to approximate the variance and confidence intervals are available.}
\item{\strong{Allowed Time Values:} Allows both continuous and integer time.}
\item{\strong{Bounded Estimates:} Estimates are guaranteed to be bounded in the 0 to 1 probability range.}
\item{\strong{Monotone Function:} Estimates are not guaranteed to be monotone.}
\item{\strong{Dependencies:} Depends on \pkg{numDeriv} and the dependencies of that package.}
}

This method is based on the proximal causal inference framework, first introduced by Miao et al. (2018) and later extended to allow for estimation of counterfactual survival curves by Ying et al. (2022). It allows the estimation of the treatment-specific counterfactual survival curve under unmeasured confounding, when the true data-generation mechanism has a particular structure. In particular, there must be an observed variable (contained in the dataset) that is a potential cause of the treatment (\code{variable}) and also unrelated to the time-to-event endpoint, except through measured confounders (\code{adjust_vars}) and a particular known but unmeasured confounder. This type of variable is called a \code{treatment_proxy}. Secondly, there must be another observed variable that directly or indirectly causes the outcome, but is unrelated to the treatment expect through measured confounders and the same known but unmeasured confounder as mentioned earlier. This type of variable is called an \code{outcome_proxy}. A better explanation is given by Zivich et al. (2023). More information on the underlying assumptions can be found in the papers listed in the references.

Ying et al. (2022) proposed two methods to utilize this kind of structure for the estimation of the counterfactual survival curve. The one implemented here relies on estimating the q-confounding bridge and the h-confounding bridge using parametric models. This essentially means that it uses both the treatment-assignment mechanism and the outcome-mechanism to adjust for confounding, similar to a regular augmented inverse probability weighted estimator.

}
\value{
Adds the following additional objects to the output of the \code{adjustedsurv} function:

\itemize{
\item \code{noncensor_cumhaz}: The estimated cumulative hazard function.
\item \code{noncensor_cumhaz_IF}: The influence function based on the estimated cumulative hazard function.
\item \code{q_bridge}: A \code{list} containing results from fitting the q-confounding bridge function.
\item \code{h_bridge}: A \code{list} containing results from fitting the h-confounding bridge function.
}
}
\references{
Andrew Ying, Yifan Cui and Eric J. Tchetgen Tchetgen (2022). "Proximal Causal Inference for Marginal Counterfactual Survival Curves". arXiv:2204.13144

Wang Miao, Zhi Geng and Eric J. Tchetgen Tchetgen (2018). "Identifying Causal Effects with Proxy Variables of an Unmeasured Confounder". In: Biometrika 105.4, pp. 987-993.

Paul N. Zivich, Stephen R. Cole, Jessie K. Edwards, Grace E. Mulholland, Bonnie E. Shook-Sa and Eric J. Tchetgen Tchetgen (2023). "Introducing Proximal Causal Inference for Epidemiologists". In: American Journal of Epidemiology 192.7, pp. 1224-1227.

Eric J. Tchetgen Tchetgen, Andrew Ying, Yifan Cui, Xu Shi and Wang Miao (2020). "An Introduction to Proximal Causal Learning". arXiv:2009.10982
}
\author{
Andrew Ying wrote all of the low-level estimation functions used to actually obtain the relevant values. Robin Denz wrote some wrapper functions around those to include this method in this package.
}

\seealso{
\code{\link[stats]{optim}}
}
\examples{
library(adjustedCurves)

#### generate some example data that fufill all assumptions ####
# code was taken from the github repository associated with the original
# paper by Ying et al. (2022): https://github.com/andrewyyp/Proximal_MSF

# simulation parameters
para_set <- list(mu_X = 1.1,
                 sigma_X = 0.75,
                 mu_U = 1.1,
                 sigma_U = 0.75,
                 alpha_A = c(0.3, 0.4, -0.6),
                 mu_Z = c(-0.2, -0.3, 0.65),
                 sigma_Z = 0.5,
                 mu_W = c(-0.6, 0.4, 0.65),
                 sigma_W = 0.5,
                 mu_T0 = c(0.1, 0.6, 0.25, 0.5),
                 mu_C = 0.2,
                 admin_C = 2
)

# small function to obtain the required data
data_gen <- function(N, para_set, a = NULL) {
  # generate X, U
  X <- para_set$mu_X + rnorm(N, 0, para_set$sigma_X)
  U <- para_set$mu_U + rnorm(N, 0, para_set$sigma_U)
  X <- pmax(X, 0)
  U <- pmax(U, 0)

  if (is.null(a)) {
    # generate A
    prop_score_0 <- 1/(1 + exp(-cbind(1, X, U) \%*\% para_set$alpha_A))
    A <- rbinom(N, 1, prop_score_0)
  } else {
    A <- rep(a, N)
  }


  # generate Z
  Z <- cbind(1, X, U) \%*\% para_set$mu_Z + rnorm(N, 0, para_set$sigma_Z)

  # generate W
  W <- cbind(1, X, U) \%*\% para_set$mu_W + rnorm(N, 0, para_set$sigma_W)


  #generate Y
  T0 <- rexp(N, rate = cbind(1, A, X, U) \%*\% para_set$mu_T0)

  C <- rexp(N, rate = para_set$mu_C)
  C <- pmin(C, para_set$admin_C)
  if (is.null(a)) {
    df <- data.frame(X, U, A, Z, W, T0 = pmin(T0, C), Delta = (T0 <= C))
  } else {
    df <- data.frame(X, U, A, Z, W, T0 = T0, Delta = rep(1, N))
  }
  return(df)
}

#### Simple example ####

set.seed(4356)
# NOTE: increase N to get more stable estimates, kept low here to pass
#       speed requirements set by CRAN
data <- data_gen(N=50, para_set=para_set)
data$A <- factor(data$A)

if (requireNamespace("numDeriv")) {

library(numDeriv)

adj <- adjustedsurv(data=data,
                    variable="A",
                    ev_time="T0",
                    event="Delta",
                    method="prox_aiptw",
                    adjust_vars="X",
                    treatment_proxy="Z",
                    outcome_proxy="W",
                    conf_int=TRUE)
plot(adj, iso_reg=TRUE)
}
}
