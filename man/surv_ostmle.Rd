\name{surv_ostmle}
\alias{surv_ostmle}

\title{
One-Step Targeted Maximum Likelihood Estimation Survival Curves
}
\description{

This page explains the details of estimating adjusted survival curves using the one-step targeted maximum likelihood methodology for single event time-to-event data (\code{method="ostmle"} in the \code{\link[=adjustedsurv]{adjustedsurv}} function). All regular arguments of the \code{adjustedsurv} function can be used. Additionally, the \code{SL.ftime}, \code{SL.ctime} and \code{SL.trt} arguments have to be specified in the \code{adjustedsurv} call. Further arguments specific to this method are listed below.

}

\arguments{
  \item{adjust_vars}{
A character vector of column names specifying variables to be used when modeling the outcome, treatment and censoring mechanism. Can be set to \code{NULL} (default), in which case all covariates are used. See details and examples.
}
  \item{SL.ftime}{
[\strong{required}] A character vector or list specification to be passed to the \code{SL.library} option in the call to \code{SuperLearner} for the outcome regression. See\code{?SuperLearner} for more information on how to specify valid \code{SuperLearner} libraries. It is expected that the wrappers used in the library will play nicely with the input variables, which will be called "trt", \code{names(adjust_vars)}, and "t".
  }
  \item{SL.ctime}{
[\strong{required}] A character vector or list specification to be passed to the \code{SL.library} argument in the call to \code{SuperLearner} for the estimate of the conditional hazard for censoring. It is expected that the wrappers used in the library will play nicely with the input variables, which will be called "trt" and \code{names(adjust_vars)}.
  }
  \item{SL.trt}{
[\strong{required}] A character vector or list specification to be passed to the \code{SL.library} argument in the call to \code{SuperLearner} for the estimate of the conditional probability of treatment. It is expected that the wrappers used in the library will play nicely with the input variables, which will be \code{names(adjust_vars)}.
  }
  \item{epsilon}{
The size of the updating step. See \code{MOSS} for more details.
}
  \item{max_num_iteration}{
The maximum number of iterations used in the updating step of the One-Step estimator.
}
  \item{psi_moss_method}{
Specifies the method used to make the pooled updating step. Can be either \code{"l1"} (Ridge Regression), \code{"l2"} (LASSO) or \code{"glm"}. For more information see documentation of \code{MOSS}.
}
  \item{tmle_tolerance}{
The tolerance used to determine whether the TMLE estimator converged or not. See \code{MOSS} for more details.
}
  \item{gtol}{
Tolerance at which estimated treatment assignment probabilities are truncated. Every propensity score bigger than 1 - \code{gtol} is set to 1 - \code{gtol} and every propensity score smaller than \code{gtol} is set to \code{gtol}. Useful when there are extreme propensity scores close to 0 or 1.
}
}
\details{

\itemize{
\item{\strong{Type of Adjustment:} Adjustments are made based on the treatment assignment mechanism, the outcome mechanism and the censoring mechanism. No models can be supplied. The adjustments are made based on \code{\link[SuperLearner]{SuperLearner}}} libraries.
\item{\strong{Doubly-Robust:} Estimates are Doubly-Robust.}
\item{\strong{Categorical groups:} Currently only two groups in \code{variable} are allowed. Must be a factor variable with exactly two levels.}
\item{\strong{Approximate Variance:} Calculations to approximate the variance and confidence intervals are available.}
\item{\strong{Allowed Time Values:} Allows only integer time.}
\item{\strong{Bounded Estimates:} Estimates are guaranteed to be bounded in the 0 to 1 probability range.}
\item{\strong{Monotone Function:} Estimates are guaranteed to be monotone.}
\item{\strong{Dependencies:} This method relies on the \pkg{survtmle}, \pkg{MASS}, \pkg{SuperLearner}, \pkg{R6} and \pkg{tidyr} packages. The code was taken directly from the \pkg{MOSS} package, which is currently not available on CRAN. Tiny changes were made but it's essentially the same. You do not have to install it but it can be installed using the following code:\cr \code{devtools::install_github("wilsoncai1992/MOSS")}. \cr If you use this method, please cite the \pkg{MOSS} package.}
}

Standard TMLE is a two-step procedure. First, initial estimates for the treatment-assignment and the outcome-mechanisms are made using loss-based learning. This is implemented here using the \code{SuperLearner} methodology. In the next step, the estimates obtained by using the outcome-mechanism model are fluctuated based on information from the treatment-assignment model. If the outcome model is already consistent, this fluctuation is very small and the estimates stay consistent. If the outcome model is biased, the fluctuation removes the bias whenever the treatment assignment model is consistent. This process is iterative and continues until a threshold is hit (either the maximum number of iterations is reached or the bias is smaller than the specified tolerance, see \code{?survtmle}).

In contrast to the estimator implemented in the \code{\link[=surv_tmle]{tmle}} method, the OSTMLE uses a LASSO or Ridge Regression model in the targeting step to fluctuate the initial estimates. Details can be found in Cai and van der Laan (2020).

As has been shown in multiple studies by Mark J. van der Laan and colleagues, this method has some desirable mathematical properties and generally performs well in appropriate scenarios. The biggest problem is however, that it is only defined for discrete (integer-valued) survival times. Simply discretizing continuous survival times only works to a certain extent and is generally discouraged.

When the sample size is large or many time points are of interest, this method will also be *very* slow. While possible to run, bootstrapping would take an enormous amount of time and is therefore discouraged.

}
\value{

Adds the following additional objects to the output of the \code{adjustedsurv} function:

\itemize{
\item{\code{psi_moss_hazard_0}: }{The iterated MOSS objects for the control group.}
\item{\code{psi_moss_hazard_1}: }{The iterated MOSS objects for the treatment group.}
}

}
\references{
Weixin Cai and Mark J. van der Laan (2020). "One-Step Targeted Maximum Likelihood Estimation for Time-To-Event Outcomes". In: Biometrics 76, pp. 722–733

David Benkeser, Marco Carone, and Peter B. Gilbert (2018). "Improved Estimation of the Cumulative Incidence of Rare Outcomes". In: Statistics in Medicine 37.2, pp. 280–293

Megan S. Schuler and Sherri Rose (2017). "Targeted Maximum Likelihood Estimation for Causal Inference in Observational Studies". In: American Journal of Epidemiology 186.1, pp. 65-73
}
\author{
The wrapper function was written by Robin Denz, the \code{MOSS} package (which this wrapper is based on) was written by Wilson Cai. See <https://github.com/wilsoncai1992/MOSS/> for more details.
}

\seealso{
\code{\link[survtmle]{survtmle}}, \code{MOSS}, \code{\link[SuperLearner]{SuperLearner}}
}
\examples{
# not run because any meaningful example is too slow
\dontrun{
library(adjustedCurves)

# simulate some data as example
sim_dat <- sim_confounded_surv(n=50, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# only works with integer time, only unbiased with small amounts of them
sim_dat$time <- round(sim_dat$time*15) + 1

# calculate adjusted survival curves, using SuperLearner but only
# using the SL.glm library. In practice you would want to use more than
# that. See ?MOSS and ?survtmle
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="ostmle",
                        adjust_vars=c("x1", "x2", "x3", "x4", "x5", "x6"),
                        SL.ftime=c("SL.glm"),
                        SL.ctim=c("SL.glm"),
                        SL.trt=c("SL.glm"))

# plot the curves
plot(adjsurv)
}
}
