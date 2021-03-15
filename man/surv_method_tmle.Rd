\name{surv_method_tmle}
\alias{surv_method_tmle}

\title{
Targeted Maximum Likelihood Estimation for Adjusted Survival Curves
}
\description{
This method utilizes the targeted maximum likelihood methodology to estimate adjusted survival curves. Similar to the AIPTW methods it posseses the doubly-robust property and can also be used to adjust for dependent censoring. It is however important to note that this method should only be used when the time variable is truly discrete and only a few distinct points in time are observed. See details for more information.
}
\usage{
surv_method_tmle(data, variable, ev_time, event, conf_int, conf_level=0.95,
                 times, adjust_vars=NULL, SL.ftime=NULL, SL.ctime=NULL,
                 SL.trt=NULL, glm.ftime=NULL, glm.ctime=NULL, glm.trt=NULL, ...)
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
A numeric vector of time points at which the survival probability should be estimated or \code{NULL}. If \code{NULL} the survival is estimated at all points in time at which an event occured. Only relevant for some methods, ignored otherwise. See details.
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
\item{\strong{Categorical groups:} Currently only two groups in \code{variable} are allowed. Must be an integer variable only containing 0 and 1.}
\item{\strong{Approximate Variance:} Calculations to approximate the variance and confidence intervals are available.}
\item{\strong{Allowed Time Values:} Allows only integer time.}
\item{\strong{Bounded Estimates:} Estimates are not guaranteed to be bounded in the 0 to 1 probability range.}
\item{\strong{Monotone Function:} Estimates are not guaranteed to be monotone.}
\item{\strong{Dependencies:} This method relies on the \code{survtmle} and \code{SuperLearner} packages.}
}

}
\value{
Returns a \code{data.frame} object containing the estimated adjusted survival probabilities for some points in time for each level of \code{variable}.
}
\references{
Megan S. Schuler and Sherri Rose. "Targeted Maximum Likelihood Estimation for Causal Inference in Observational Studies". In: American Journal of Epidemiology 186.1 (2017), pp. 65-73.


}
\author{
The wrapper function was written by Robin Denz, the \code{survtmle} package (which this wrapper is based on) was written by David Benkeser and Nima Hejazi. See \code{?survtmle} for more details.
}

\seealso{
\code{\link{survtmle}}, \code{\link[SuperLearner]{SuperLearner}}, \code{\link{glm}}
}
\examples{
\dontrun{
# simulate some data as example (needs a binary integer "variable")
sim_dat <- sim_confounded_surv(n=500, max_t=1.2)

# only works with integer time, only unbiased with small amounts of them
sim_dat$time <- round(sim_dat$time*15) + 1

# calculate adjusted survival curves, using SuperLearner but only
# using the SL.glm library. In practice you would want to use more than
# that. See ?survtmle
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="tmle",
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
