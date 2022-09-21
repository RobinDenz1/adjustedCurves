\name{cif_tmle}
\alias{cif_tmle}

\title{
Targeted Maximum Likelihood Estimation for CIFs
}
\description{

This page explains the details of estimating adjusted cumulative incidence functions using the targeted maximum likelihood methodology in a competing risks setting (\code{method="tmle"} in the \code{\link[=adjustedcif]{adjustedcif}} function). All regular arguments of the \code{adjustedcif} function can be used. Additionally, you have to specify either \code{SL.ftime} or \code{glm.ftime}, \code{SL.ctime} or \code{glm.ctime} and \code{SL.trt} or \code{glm.trt} in the \code{adjustedcif} call. Further arguments specific to this method are listed below.

}

\arguments{
  \item{adjust_vars}{
A character vector of column names specifying variables to be used when modeling the outcome, treatment and censoring mechanism. Can be set to \code{NULL} (default), in which case all covariates are used. See details and examples.
  }
  \item{SL.ftime}{
A character vector or list specification to be passed to the \code{SL.library} option in the call to \code{SuperLearner} for the outcome regression. See\code{?SuperLearner} for more information on how to specify valid \code{SuperLearner} libraries. It is expected that the wrappers used in the library will play nicely with the input variables, which will be called "trt", \code{names(adjust_vars)}, and "t".
  }
  \item{SL.ctime}{
A character vector or list specification to be passed to the \code{SL.library} argument in the call to \code{SuperLearner} for the estimate of the conditional hazard for censoring. It is expected that the wrappers used in the library will play nicely with the input variables, which will be called "trt" and \code{names(adjust_vars)}.
  }
  \item{SL.trt}{
A character vector or list specification to be passed to the \code{SL.library} argument in the call to \code{SuperLearner} for the estimate of the conditional probability of treatment. It is expected that the wrappers used in the library will play nicely with the input variables, which will be \code{names(adjust_vars)}.
  }
  \item{glm.ftime}{
A character specification of the right-hand side of the equation passed to the \code{formula} option of a call to \code{glm} for the outcome regression. Ignored if \code{SL.ftime} is not equal to \code{NULL}. Use "trt" to specify the treatment in this formula (see examples). The formula can additionally include any variables found in \code{names(adjust_vars)}.
  }
  \item{glm.ctime}{
A character specification of the right-hand side of the equation passed to the \code{formula} option of a call to \code{glm} for the estimate of the conditional hazard for censoring.  Ignored if \code{SL.ctime} is not equal to \code{NULL}. Use "trt" to specify the treatment in this formula (see examples).  The formula can additionally include any variables found in \code{names(adjust_vars)}.
  }
  \item{glm.trt}{
A character specification of the right-hand side of the equation passed to the \code{formula} option of a call to \code{glm} for the estimate of the conditional probability of treatment. Ignored if \code{SL.trt} is not equal to \code{NULL}. The formula can include any variables found in \code{names(adjust_vars)}.
  }
  \item{...}{
Additional arguments passed to \code{survtmle}.
  }
}
\details{

\itemize{
\item{\strong{Type of Adjustment:} Adjustments are made based on the treatment assignment mechanism, the outcome mechanism and the censoring mechanism. No models can be supplied. The adjustments are made based on \code{\link[SuperLearner]{SuperLearner}}} libraries or using the \code{glm} arguments.
\item{\strong{Doubly-Robust:} Estimates are Doubly-Robust.}
\item{\strong{Categorical groups:} Currently only two groups in \code{variable} are allowed. Must be a factor variable with exactly two levels.}
\item{\strong{Approximate Variance:} Calculations to approximate the variance and confidence intervals are available.}
\item{\strong{Allowed Time Values:} Allows only integer time.}
\item{\strong{Bounded Estimates:} Estimates are guaranteed to be bounded in the 0 to 1 probability range.}
\item{\strong{Monotone Function:} Estimates are not guaranteed to be monotone.}
\item{\strong{Dependencies:} This method relies on the \pkg{survtmle} and \pkg{SuperLearner} packages.}
}

TMLE is a two-step procedure. First, initial estimates for the treatment-assignment and the outcome-mechanisms are made using loss-based learning. This is implemented here using the \code{SuperLearner} methodology. In the next step, the estimates obtained by using the outcome-mechanism model are fluctuated based on information from the treatment-assignment model. If the outcome model is already consistent, this fluctuation is very small and the estimates stay consistent. If the outcome model is biased, the fluctuation removes the bias whenever the treatment assignment model is consistent. This process is iterative and continues until a threshold is hit (either the maximum number of iterations is reached or the bias is smaller than the specified tolerance, see \code{?survtmle}).

As has been shown in multiple studies by Mark J. van der Laan and colleagues, this method has some desirable mathematical properties and generally performs well in appropriate scenarios. The biggest problem is however, that it is only defined for discrete (integer-valued) survival times. Simply discretizing continuous survival times only works to a certain extent and is generally discouraged.

When the sample size is large or many time points are of interest, this method will also be *very* slow. While possible to run, bootstrapping would take an enormous amount of time and is therefore discouraged.

}
\value{
Adds the following additional objects to the output of the \code{adjustedcif} function:

\itemize{
\item{\code{survtmle_object}: }{The object created using the \code{\link[survtmle]{survtmle}} function.}
\item{\code{survtmle.timepoints_object}: }{The object created using the \code{survtmle.timepoints} function.}
}
}
\references{
Megan S. Schuler and Sherri Rose (2017). "Targeted Maximum Likelihood Estimation for Causal Inference in Observational Studies". In: American Journal of Epidemiology 186.1, pp. 65-73

David Benkeser, Marco Carone, and Peter B. Gilbert (2018). "Improved Estimation of the Cumulative Incidence of Rare Outcomes". In: Statistics in Medicine 37.2, pp. 280-293
}
\author{
The wrapper function was written by Robin Denz, the \pkg{survtmle} package (which this wrapper is based on) was written by David Benkeser and Nima Hejazi. See \code{?survtmle} for more details.
}

\seealso{
\code{\link[survtmle]{survtmle}}, \code{\link[SuperLearner]{SuperLearner}}, \code{\link{glm}}
}
\examples{
# not run because any meaningful example is too slow
\donttest{
library(adjustedCurves)
library(survtmle)

set.seed(42)

# simulate some data as example
sim_dat <- sim_confounded_crisk(n=30, max_t=5)
sim_dat$group <- as.factor(sim_dat$group)

# only works with integer time, only unbiased with small amounts of them
sim_dat$time <- round(sim_dat$time*15) + 1

# calculate adjusted CIFs, using SuperLearner but only
# using the SL.glm library. In practice you would want to use more than
# that. See ?survtmle
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      cause=1,
                      method="tmle",
                      adjust_vars=c("x1", "x2", "x3", "x4", "x5", "x6"),
                      SL.ftime=c("SL.glm"),
                      SL.ctim=c("SL.glm"),
                      SL.trt=c("SL.glm"))

# plot the curves
plot(adjcif, iso_reg=TRUE)
}
}
