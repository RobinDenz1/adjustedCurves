\name{surv_ostmle}
\alias{surv_ostmle}

\title{
One-Step Targeted Maximum Likelihood Estimation Survival Curves
}
\description{
This method utilizes the one-step targeted maximum likelihood methodology to estimate adjusted survival curves. It is very similar to the \code{\link[=surv_tmle]{tmle}} method, with some difference in the estimation techniques. Similar to the AIPTW methods it posseses the doubly-robust property and can also be used to adjust for dependent censoring. It is however important to note that this method should only be used when the time variable is truly discrete and only a few distinct points in time are observed. See details for more information.
}
\usage{
surv_ostmle(data, variable, ev_time, event, conf_int,
            conf_level=0.95, times, adjust_vars=NULL,
            SL.ftime=NULL, SL.ctime=NULL, SL.trt=NULL,
            epsilon=1, max_num_iteration=100,
            psi_moss_method="l2", tmle_tolerance=NULL,
            gtol=1e-3)
}

\arguments{
  \item{data}{
A \code{data.frame} object containing the needed time-to-event data in standard format.
}
  \item{variable}{
A character string specifying the variable by which the survival curves should be grouped. Must be a valid columname of \code{data}.
}
  \item{ev_time}{
A character string specifying the variable indicating the time-to-event or time-to-censoring. Must be a valid columname of \code{data}.
}
  \item{event}{
A character string specifying the binary event indicator. Must be a valid columname of \code{data}.
}
  \item{conf_int}{
A logical variable, indicating whether the asymptotic variances and confidence intervals of the survival probabilities should be calculated.
  }
  \item{conf_level}{
A number specifying the confidence level of asymptotic confidence intervals.
  }
  \item{times}{
A numeric vector of time points at which the survival probability should be estimated.
}
  \item{adjust_vars}{
A character vector of column names specifying variables to be used when modeling the outcome, treatment and censoring mechanism. See details and examples.
}
  \item{SL.ftime}{
A character vector or list specification to be passed to the \code{SL.library} option in the call to \code{SuperLearner} for the outcome regression. See\code{?SuperLearner} for more information on how to specify valid \code{SuperLearner} libraries. It is expected that the wrappers used in the library will play nicely with the input variables, which will be called "trt", \code{names(adjust_vars)}, and "t".
  }
  \item{SL.ctime}{
A character vector or list specification to be passed to the \code{SL.library} argument in the call to \code{SuperLearner} for the estimate of the conditional hazard for censoring. It is expected that the wrappers used in the library will play nicely with the input variables, which will be called "trt" and \code{names(adjust_vars)}.
  }
  \item{SL.trt}{
A character vector or list specification to be passed to the \code{SL.library} argument in the call to \code{SuperLearner} for the estimate of the conditional probability oftreatment. It is expected that the wrappers used in the library will play nicely with the input variables, which will be \code{names(adjust_vars)}.
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
Tolerance at which estimated treatment assignment probabilities are truncated.
}
}
\details{

\itemize{
\item{\strong{Type of Adjustment:} Adjustments are made based on the treatment assignment mechanism, the outcome mechanism and the censoring mechanism. No models can be supplied. The adjustments are made based on \code{\link[SuperLearner]{SuperLearner}}} libraries.
\item{\strong{Doubly-Robust:} Estimates are Doubly-Robust.}
\item{\strong{Categorical groups:} Currently only two groups in \code{variable} are allowed. Must be an integer variable only containing 0 and 1.}
\item{\strong{Approximate Variance:} Calculations to approximate the variance and confidence intervals are available.}
\item{\strong{Allowed Time Values:} Allows only integer time.}
\item{\strong{Bounded Estimates:} Estimates are guaranteed to be bounded in the 0 to 1 probability range.}
\item{\strong{Monotone Function:} Estimates are guaranteed to be monotone.}
\item{\strong{Dependencies:} This method relies on the \code{survtmle}, \code{MASS}, \code{SuperLearner}, \code{R6} and \code{tidyr} packages. The code was taken directly from the \code{MOSS} package, which is currently not available on CRAN. We made tiny changes to it but it's essentially the same. You do not have to install it but it can be installed using the following code: \code{devtools::install_github("wilsoncai1992/MOSS")}. If you use this method, please cite the \code{MOSS} package.}
}

Standard TMLE is a two-step procedure. First, initial estimates for the treatment-assignment and the outcome-mechanisms are made using loss-based learning. This is implemented here using the \code{SuperLearner} methodology. In the next step, the estimates obtained by using the outcome-mechanism model are fluctuated based on information from the treatment-assignment model. If the outcome model is already consistent, this fluctuation is very small and the estimates stay consistent. If the outcome model is biased, the fluctuation removes the bias whenever the treatment assignment model is consistent. This process is iterative and continues until a threshhold is hit (either the maximum number of iterations is reached or the bias is smaller than the specified tolerance, see \code{?survtmle}).

In contrast to the estimator implemented in the \code{\link[=surv_tmle]{tmle}} method, the OSTMLE uses a LASSO or Ridge Regression model in the targeting step to fluctuate the initial estimates. Details can be found in Cai and van der Laan (2020).

As has been shown in multiple studies by Mark J. van der Laan and colleagues, this method has some desireable mathematical properties and generally performs well in appropriate scenarios. The biggest problem is however, that it is only defined for discrete (integer-valued) survival times. Simply discretizing continuous survival times only works to a certain extent and is generally discouraged. If TMLE estimates for continuous survival data are the goal, the user should use the \code{\link[=surv_tmle_pseudo]{tmle_pseudo}} method instead.

When the sample size is large or many time points are of interest, this method will also be *very* slow. While possible to run, bootstrapping would take an enormous amount of time and is therefore discouraged.

}
\value{
Returns a \code{data.frame} object containing the estimated adjusted survival probabilities for some points in time for each level of \code{variable}.
}
\references{
Weixin Cai and Mark J. van der Laan. "One-Step Targeted Maximum Likelihood Estimation for Time-To-Event Outcomes". In: Biometrics 76 (2020), pp. 722-733.
}
\author{
The wrapper function was written by Robin Denz, the \code{MOSS} package (which this wrapper is based on) was written by Wilson Cai. See <https://github.com/wilsoncai1992/MOSS/> for more details.
}

\seealso{
\code{\link[survtmle]{survtmle}}, \code{MOSS}, \code{\link[SuperLearner]{SuperLearner}}
}
\examples{
\dontrun{
# simulate some data as example (needs a binary integer "variable")
sim_dat <- sim_confounded_surv(n=500, max_t=1.2)

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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
