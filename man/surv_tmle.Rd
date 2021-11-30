\name{surv_tmle}
\alias{surv_tmle}

\title{
Targeted Maximum Likelihood Estimation for Adjusted Survival Curves
}
\description{

This page explains the details of estimating adjusted survival curves using the targeted maximum likelihood methodology for single event time-to-event data. The function is called internally when using \code{method="tmle"} in the \code{\link[=adjustedsurv]{adjustedsurv}} function. It is highly recommended to use \code{\link[=adjustedsurv]{adjustedsurv}} instead of calling this function directly.

}
\usage{
surv_tmle(data, variable, ev_time, event, conf_int,
          conf_level=0.95, times, adjust_vars=NULL,
          SL.ftime=NULL, SL.ctime=NULL,
          SL.trt=NULL, glm.ftime=NULL,
          glm.ctime=NULL, glm.trt=NULL, ...)
}

\arguments{
  \item{data}{
A \code{data.frame} object containing the needed time-to-event data in standard format.
  }
  \item{variable}{
A character string specifying the variable by which the survival curves should be grouped. Must be a valid column name of \code{data}. Can be either a factor variable or a binary integer variable. If a factor variable is supplied, the first level is assumed to correspond to the control group. Either way only two treatments are allowed.
  }
  \item{ev_time}{
A character string specifying the variable indicating the time-to-event or time-to-censoring. Must be a valid column name of \code{data}.
  }
  \item{event}{
A character string specifying the binary event indicator. Must be a valid column name of \code{data}.
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
Additional arguments passed to \code{\link[survtmle]{survtmle}}.
  }
}
\details{

\itemize{
\item{\strong{Type of Adjustment:} Adjustments are made based on the treatment assignment mechanism, the outcome mechanism and the censoring mechanism. No models can be supplied. The adjustments are made based on \code{\link[SuperLearner]{SuperLearner}}} libraries or using the \code{glm} arguments.
\item{\strong{Doubly-Robust:} Estimates are Doubly-Robust.}
\item{\strong{Categorical groups:} Currently only two groups in \code{variable} are allowed. Must be an integer variable only containing 0 and 1 or a factor variable with exactly two levels.}
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
Returns a \code{list} object containing a \code{data.frame} with the estimated adjusted survival probabilities for some points in time for each level of \code{variable}, the object created using the \code{\link[survtmle]{survtmle}} function and the object created using the \code{survtmle.timepoints} function.
}
\references{
Ori M. Stitelman and Mark J. van der Laan (2010). "Collaborative Targeted Maximum Likelihood for Time to Event Data". In: The International Journal of Biostatistics 6.1

David Benkeser, Marco Carone, and Peter B. Gilbert (2018). "Improved Estimation of the Cumulative Incidence of Rare Outcomes". In: Statistics in Medicine 37.2, pp. 280-293

Megan S. Schuler and Sherri Rose (2017). "Targeted Maximum Likelihood Estimation for Causal Inference in Observational Studies". In: American Journal of Epidemiology 186.1, pp. 65-73
}
\author{
The wrapper function was written by Robin Denz, the \code{survtmle} package (which this wrapper is based on) was written by David Benkeser and Nima Hejazi. See \code{?survtmle} for more details.
}

\seealso{
\code{\link[survtmle]{survtmle}}, \code{\link[SuperLearner]{SuperLearner}}, \code{\link{glm}}
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
