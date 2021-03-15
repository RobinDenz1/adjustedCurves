\name{cif_method_matching}
\alias{cif_method_matching}

\title{
Using Matching to Calculate Adjusted CIFs
}
\description{
This method works by first matching controls and treated patients using a previously estimated propensity score. Using this new matched population, a simple Aalen-Johansen estimator is used to calculate the adjusted CIFs.
}
\usage{
cif_method_matching(data, variable, ev_time, event, cause, conf_int,
                    conf_level=0.95, treatment_model, ...)
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
A logical variable, indicating whether the asymptotic variances and confidence intervals of the CIFs should be calculated.
  }
  \item{conf_level}{
A number specifying the confidence level of asymptotic and/or bootstrap confidence intervals.
  }
  \item{treatment_model}{
Must be either a model object with \code{variable} as response variable or a vector of previously estimated propensity scores.
}
  \item{...}{
Further arguments passed to the \code{Match} function of the \code{Matching} Package.
  }
}
\details{

\itemize{
\item{\strong{Type of Adjustment:} Requires a model describing the treatment assignment mechanism. This must be either a \code{\link{glm}} object or a vector of propensity scores.}
\item{\strong{Doubly-Robust:} Estimates are not Doubly-Robust.}
\item{\strong{Categorical groups:} Only two groups in \code{variable} are allowed. Must be an integer variable containing only 0 and 1.}
\item{\strong{Approximate Variance:} Calculations to approximate the variance and confidence intervals are currently not available. Bootstrapping is also discouraged.}
\item{\strong{Allowed Time Values:} Allows both continuous and integer time.}
\item{\strong{Bounded Estimates:} Estimates are guaranteed to be bounded in the 0 to 1 probability range.}
\item{\strong{Monotone Function:} Estimates are guaranteed to be monotone.}
\item{\strong{Dependencies:} This method relies on both the \code{Matching} and the \code{cmprsk} packages.}
}

WARUM KEIN BOOTSTRAP?

}
\value{
Returns a \code{data.frame} object containing the estimated adjusted CIFs for \code{cause} at some points in time for each level of \code{variable}.
}
\references{
Angela Winnett and Peter Sasieni. "Adjusted Nelson-Aalen Estimates with Retrospective Matching". In: Journal of the American Statistical Association 97.457 (2002), pp. 245-256.

Peter C. Austin, Neal Thomas, and Donald B. Rubin. "Covariate-Adjusted Survival Analyses in Propensity-Score Matched Samples: Imputing Potential Time-To-Event Outcomes". In: Statistical Methods in Medical Research 29.3 (2020), pp. 728-751.
}
\author{
The wrapper function was written by Robin Denz, the \code{Match} function (which this wrapper is build around) was written by other people. See \code{?Match} for more details.
}

\seealso{
\code{\link[Matching]{Match}}, \code{\link[cmprsk]{cuminc}}
}
\examples{
\dontrun{
# simulate some data as example (needs a binary integer "variable")
sim_dat <- sim_confounded_surv(n=500, max_t=1.2)

# estimate treatment assignment model
glm_mod <- glm(group ~ x1 + x2 + x4 + x6, data=sim_dat, family="binomial")

# calculate adjusted survival curves
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="matching",
                        treatment_model=glm_mod)

# Alternatively, supply the propensity score directly
# Here we use the logistic regression to calculate it, so we get
# exactly the same result. The propensity score can be calculated in
# any other way in practice, allowing flexibility
ps_score <- glm_mod$fitted.values

adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="matching",
                        treatment_model=ps_score)

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
