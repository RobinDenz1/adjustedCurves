\name{surv_km}
\alias{surv_km}

\title{
Group-Specific Kaplan-Meier Survival Curves
}
\description{

This page explains the details of estimating group-specific Kaplan-Meier curves for single event time-to-event data. The function is called internally when using \code{method="km"} in the \code{\link[=adjustedsurv]{adjustedsurv}} function.

Calculates standard Kaplan-Meier survival curves, stratified by the group variable.
NO adjustment for any confounders is made. This function is included only for reference and should not be used when confounder adjusted survival curves are desired.

}
\usage{
surv_km(data, variable, ev_time, event, conf_int, conf_level=0.95,
        times=NULL, conf_type="log")
}

\arguments{
  \item{data}{
A \code{data.frame} object containing the needed time-to-event data in standard format.
}
  \item{variable}{
A character string specifying the variable by which the survival curves should be grouped. Must be a valid columname of \code{data}. The variable specified should be a factor variable.
}
  \item{ev_time}{
A character string specifying the variable indicating the time-to-event or time-to-censoring. Must be a valid columname of \code{data}.
}
  \item{event}{
A character string specifying the binary event indicator. Must be a valid columname of \code{data}.
}
  \item{conf_int}{
A logical variable, indicating whether the asymptotic variances and confidence intervals of the survival probabilities should be calculated. For the available types of confidence intervals see \code{?survfit}.
  }
  \item{conf_level}{
A number specifying the confidence level of asymptotic and/or bootstrap confidence intervals.
  }
  \item{times}{
A numeric vector of time points for which to return the survival probability estimates or \code{NULL} (default). If \code{NULL} the whole curve will be estimated according to the usual Kaplan-Meier procedure. When custom time points are specified the whole curve is estimated anyways and the time specific estimates are simply read off that curve.
}
  \item{conf_type}{
The type of confidence interval that should be calculated. Has to be a character string, passed to the \code{conf.type} argument in the \code{survfit} function. Defaults to \code{"log"}, which is also the default in \code{survfit}.
  }
}
\details{

\itemize{
\item{\strong{Type of Adjustment:} NO adjustments are made. This is just a stratified Kaplan-Meier estimator.}
\item{\strong{Doubly-Robust:} Estimates are not Doubly-Robust.}
\item{\strong{Categorical groups:} Any number of levels in \code{variable} are allowed. Must be a factor variable.}
\item{\strong{Approximate Variance:} Calculations to approximate the variance and confidence intervals are available.}
\item{\strong{Allowed Time Values:} Allows both continuous and integer time.}
\item{\strong{Bounded Estimates:} Estimates are guaranteed to be bounded in the 0 to 1 probability range.}
\item{\strong{Monotone Function:} Estimates are guaranteed to be monotone.}
\item{\strong{Dependencies:} This method relies on the the \code{\link{survival}} package.}
}

}
\value{
Returns a \code{list} object containing a \code{data.frame} with the estimated adjusted survival probabilities for some points in time for each level of \code{variable} and the original \code{survfit} object.
}
\references{
E. L. Kaplan and Paul Meier. "Nonparametric Estimation from Incomplete Observations". In: Journal of the American Statistical Association 53.282 (1958), pp. 457-481.
}
\author{
The wrapper function was written by Robin Denz, the \code{survfit} function (which this wrapper is build around) was written by other people. See \code{?survfit} for more details.
}

\seealso{
\code{\link{survival}}, \code{\link{survfit}}
}
\examples{
\dontrun{
library(survival)

# simulate some data as example
sim_dat <- sim_confounded_surv(n=500, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# calculate un-adjusted kaplan-meier survival curves
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="km")

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
