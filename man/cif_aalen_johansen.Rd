\name{cif_aalen_johansen}
\alias{cif_aalen_johansen}

\title{
Group-Specific Aalen-Johansen CIFs
}
\description{
Calculates standard Aalen-Johansen cumulative incidence functions, stratified by the group variable.
No adjustment for any confounders is made. This function is included only for reference and should not be used when confounder adjusted CIFs are desired.
}
\usage{
cif_aalen_johansen(data, variable, ev_time, event, cause, conf_int,
                   conf_level=0.95, ...)
}

\arguments{
  \item{data}{
A \code{data.frame} object containing the needed time-to-event data in standard format.
}
  \item{variable}{
A character string specifying the variable by which the CIFs should be grouped. Must be a valid columname of \code{data}.
}
  \item{ev_time}{
A character string specifying the variable indicating the time-to-event or time-to-censoring. Must be a valid columname of \code{data}.
}
  \item{event}{
A character string specifying the numeric event indicator. The censoring indicator should be coded as 0 and all other events of interest as 1, 2, etc. Must be a valid columname of \code{data}.
  }
  \item{cause}{
The cause of interest for which the cumulative incidence functions should be calculated. Should be a number that appears in the \code{event} column of \code{data}.
  }
  \item{conf_int}{
A logical variable, indicating whether the asymptotic variances and confidence intervals of the CIFs should be calculated. For the available types of confidence intervals see \code{?cuminc}.
  }
  \item{conf_level}{
A number specifying the confidence level of asymptotic confidence intervals.
  }
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
\item{\strong{Dependencies:} This method relies on the the \code{cmprsk} package.}
}

This function is just a convenient wrapper around the \code{\link[cmprsk]{cuminc}} function. See \code{?cuminc} or the cited literature for more details.

}
\value{
Returns a \code{data.frame} object containing the estimated CIFs for some points in time for each level of \code{variable}.
}
\references{
Odd O. Aalen and Søren Johansen. "An Empirical Transition Matrix for Non-Homogeneous Markov Chains Based on Censored Observations". In: Scandinavian Journal of Statistics 5.3 (1978), pp. 141–150.
}
\author{
The wrapper function was written by Robin Denz, the \code{cuminc} function (which this wrapper is build around) was written by other people. See \code{?cuminc} for more details.
}

\seealso{
\code{\link[cmprsk]{cuminc}}
}
\examples{
\dontrun{
library(cmprsk)

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
