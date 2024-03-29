\name{cif_aalen_johansen}
\alias{cif_aalen_johansen}

\title{
Group-Specific Aalen-Johansen CIFs
}
\description{

This page explains the details of estimating standard Aalen-Johansen cumulative incidence functions, stratified by the group variable (\code{method="aalen_johansen"} in the \code{\link[=adjustedcif]{adjustedcif}} function). All regular arguments of the \code{adjustedcif} function can be used. Further arguments specific to this method are listed below.

NO adjustment for any confounders are made. This function is included only for reference and should not be used when confounder adjusted CIFs are desired.
}

\arguments{
  \item{...}{
Further arguments passed to \code{\link[cmprsk]{cuminc}}.
  }
}
\details{

\itemize{
\item{\strong{Type of Adjustment:} NO adjustments are made. This is just a stratified Aalen-Johansen estimator.}
\item{\strong{Doubly-Robust:} Estimates are not Doubly-Robust.}
\item{\strong{Categorical groups:} Any number of levels in \code{variable} are allowed. Must be a factor variable.}
\item{\strong{Approximate Variance:} Calculations to approximate the variance and confidence intervals are available.}
\item{\strong{Allowed Time Values:} Allows both continuous and integer time.}
\item{\strong{Bounded Estimates:} Estimates are guaranteed to be bounded in the 0 to 1 probability range.}
\item{\strong{Monotone Function:} Estimates are guaranteed to be monotone.}
\item{\strong{Dependencies:} This method relies on the the \pkg{cmprsk} package.}
}

This function is just a convenient wrapper around the \code{\link[cmprsk]{cuminc}} function. See \code{?cuminc} or the cited literature for more details.

}
\value{
Adds the following additional objects to the output of the \code{adjustedsurv} function:

\itemize{
\item \code{cuminc_object}: The object returned by the \code{cuminc} function.
}
}
\references{
Odd O. Aalen and Søren Johansen (1978). "An Empirical Transition Matrix for Non-Homogeneous Markov Chains Based on Censored Observations". In: Scandinavian Journal of Statistics 5.3, pp. 141-150
}
\author{
The wrapper function was written by Robin Denz, the \code{cuminc} function (which this wrapper is build around) was written by other people. See \code{?cuminc} for more details.
}

\seealso{
\code{\link[cmprsk]{cuminc}}
}
\examples{
library(adjustedCurves)

if (requireNamespace("cmprsk")) {

library(cmprsk)

set.seed(42)

# simulate some data as example
sim_dat <- sim_confounded_crisk(n=50, max_t=5)
sim_dat$group <- as.factor(sim_dat$group)

# calculate un-adjusted aalen-johansen estimates
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      cause=1,
                      method="aalen_johansen")

# plot the curves
plot(adjcif)
}
}
