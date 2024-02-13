\name{surv_km}
\alias{surv_km}

\title{
Group-Specific Kaplan-Meier Survival Curves
}
\description{

This page explains the details of estimating group-specific Kaplan-Meier curves for single event time-to-event data (\code{method="km"} in the \code{\link[=adjustedsurv]{adjustedsurv}} function). All regular arguments of the \code{adjustedsurv} function can be used. Further arguments specific to this method are listed below.

Calculates standard Kaplan-Meier survival curves, stratified by the group variable.
NO adjustment for any confounders is made. This function is included only for reference and should not be used when confounder adjusted survival curves are desired.

}

\arguments{
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
\item{\strong{Dependencies:} This method relies on the the \pkg{survival} package.}
}

}
\value{
Adds the following additional objects to the output of the \code{adjustedsurv} function:

\itemize{
\item \code{survfit_object}: The \code{survfit} object used to calculate the Kaplan-Meier curves.
}
}
\references{
E. L. Kaplan and Paul Meier (1958). "Nonparametric Estimation from Incomplete Observations". In: Journal of the American Statistical Association 53.282, pp. 457-481
}
\author{
The wrapper function was written by Robin Denz, the \code{survfit} function (which this wrapper is build around) was written by other people. See \code{?survfit} for more details.
}

\seealso{
\code{\link{survfit}}
}
\examples{
library(adjustedCurves)
library(survival)

set.seed(42)

# simulate some data as example
sim_dat <- sim_confounded_surv(n=50, max_t=1.2)
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
