\name{surv_strat_gregory}
\alias{surv_strat_gregory}

\title{
Adjusted Survival Curves for Categorical Confounders using the Method by Gregory (1988)
}
\description{

This page explains the details of estimating confounder-adjusted survival curves using a weighted average of stratified Kaplan-Meier estimates using the method described in Gregory (1988). The function is called internally when using \code{method="strat_gregory"} in the \code{\link[=adjustedsurv]{adjustedsurv}} function. It is highly recommended to use \code{\link[=adjustedsurv]{adjustedsurv}} instead of calling this function directly.

}
\usage{
surv_strat_gregory(data, variable, ev_time, event,
                   conf_int=FALSE, conf_level=0.95,
                   times=NULL, adjust_vars, na.rm=FALSE)
}

\arguments{
  \item{data}{
A \code{data.frame} object containing the needed time-to-event data in standard format.
}
  \item{variable}{
A character string specifying the variable by which the survival curves should be grouped. Must be a valid column name of \code{data}. The variable specified should be a factor variable.
}
  \item{ev_time}{
A character string specifying the variable indicating the time-to-event or time-to-censoring. Must be a valid column name of \code{data}.
}
  \item{event}{
A character string specifying the binary event indicator. Must be a valid column name of \code{data}.
}
  \item{conf_int}{
Currently not used because asymptotic variance calculations for this estimator have not been implemented.
  }
  \item{conf_level}{
A number specifying the confidence level of asymptotic confidence intervals.
  }
  \item{times}{
A numeric vector of time points at which the survival probability should be estimated or \code{NULL} (default). If \code{NULL}, the full survival curves are estimated using the standard method. If \code{times} is numeric, the full curves are still estimated and the respective values are simply read off the curves.
}
  \item{adjust_vars}{
A single string or character vector specifying column names in data for which the survival curves should be adjusted for. The variables specified can be integers, factors or characters. Only categorical variables can be used with this method. See details.
  }
  \item{na.rm}{
Whether to remove the rows from the output at which no survival probability could be estimated. The reason for the occurence of \code{NA} values is given in the details.
  }
}
\details{

\itemize{
\item{\strong{Type of Adjustment:} The survival curves are adjusted by taking a weighted average of stratified Kaplan-Meier estimates. This only works for categorical confounders. See below for more information.}
\item{\strong{Doubly-Robust:} Estimates are not Doubly-Robust.}
\item{\strong{Categorical groups:} Any number of levels in \code{variable} are allowed. Must be a factor variable.}
\item{\strong{Approximate Variance:} Calculations to approximate the variance and confidence intervals are not available.}
\item{\strong{Allowed Time Values:} Allows both continuous and integer time.}
\item{\strong{Bounded Estimates:} Estimates are guaranteed to be bounded in the 0 to 1 probability range.}
\item{\strong{Monotone Function:} Estimates are guaranteed to be monotone.}
\item{\strong{Dependencies:} This method has no additional dependencies.}
}

This is one of the older adjustment methods desribed in the literature. It only works for categorical confounders. If adjustments for continuous confounders are desired, the user needs to explicitly categorize the continuous confounders. It is recommended to use one of the other methods implemented in this package in that case. The method works exactly as desribed in Gregory (1988). TODO. A more detailed description can be found in the original article.

If a character vector is supplied in the \code{adjust_vars} argument, the Kaplan-Meier estimates are created for each combination of all supplied variables. If the sample size is small and/or there are many levels in these variables, the estimates can become unstable or undefined. Because it is a weighted average of Kaplan-Meier curves, estimates for this method are only defined for points in time with a valied Kaplan-Meier estimate in all strata. Continuing the example from above, if the Kaplan-Meier curve of the strata "Treatment + male" only extends to t = 100, it will be impossible to estimate the adjusted survival curve for t > 100 using this method.

}
\value{
Returns a \code{list} object containing a \code{data.frame} with the estimated adjusted survival probabilities for some points in time for each level of \code{variable}.
}
\references{
Gregory, W. M. Adjusting Survival Curves for Imbalances in Prognostic Factors. British Journal of Cancer, 1988, 58, 202-204
}
\author{
Robin Denz
}

\seealso{
\code{\link[adjustedCurves]{adjustedsurv}}
}
\examples{
\dontrun{
library(survival)

# simulate some data as example
sim_dat <- sim_confounded_surv(n=500, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# adjust survival curves for some categorical confounders
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="strat_gregory",
                        adjust_vars=c("x1", "x3"),
                        conf_int=FALSE)

# plot the curves
plot(adjsurv)
}
}
