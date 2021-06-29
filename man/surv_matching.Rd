\name{surv_matching}
\alias{surv_matching}

\title{
Using Matching to Calculate Adjusted Survival Curves
}
\description{

This page explains the details of estimating adjusted survival curves using propensity-score matching for single event time-to-event data. The function is called internally when using \code{method="matching"} in the \code{\link[=adjustedsurv]{adjustedsurv}} function.

This method works by first matching controls and treated patients using a previously estimated propensity score. Using this new matched population, a simple Kaplan-Meier estimator is used to calculate the adjusted survival curves. Simulation studies have shown that this method, while unbiased if used appropriately, is significantly less efficient than other methods.

}
\usage{
surv_matching(data, variable, ev_time, event, conf_int,
              conf_level=0.95, times, treatment_model, stabilize=T,
              gtol=0.001, ...)
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
A numeric vector of time points for which to return the survival probability estimates or \code{NULL} (default). If \code{NULL} the whole curve will be estimated according to the usual Kaplan-Meier procedure. When custom time points are specified the whole curve is estimated anyways and the time specific estimates are simply read off that curve.
}
  \item{treatment_model}{
Must be either a model object with \code{variable} as response variable or a vector of previously estimated propensity scores.
}
  \item{stabilize}{
Whether to stabilize the weights or not. Is set to \code{TRUE} by default. Stabilizing has no effect on the estimates, only on the asymptotic variance calculations and confidence intervals. Without stabilizing, those estimates are biased downward.
  }
  \item{gtol}{
Tolerance at which estimated treatment assignment probabilities are truncated. Every propensity score bigger than 1 - \code{gtol} is set to 1 - \code{gtol} and every propensity score smaller than \code{gtol} is set to \code{gtol}. Useful when there are extreme propensity scores close to 0 or 1.
  }
  \item{...}{
Further arguments passed to the \code{Match} function of the \code{Matching} Package.
  }
}
\details{

\itemize{
\item{\strong{Type of Adjustment:} Requires a model describing the treatment assignment mechanism. This must be either a \code{\link{glm}} object.}
\item{\strong{Doubly-Robust:} Estimates are not Doubly-Robust.}
\item{\strong{Categorical groups:} Only two groups in \code{variable} are allowed. Must be an integer variable containing only 0 and 1.}
\item{\strong{Approximate Variance:} Calculations to approximate the variance and confidence intervals are currently not available. Bootstrapping is also discouraged.}
\item{\strong{Allowed Time Values:} Allows both continuous and integer time.}
\item{\strong{Bounded Estimates:} Estimates are guaranteed to be bounded in the 0 to 1 probability range.}
\item{\strong{Monotone Function:} Estimates are guaranteed to be monotone.}
\item{\strong{Dependencies:} This method relies on the \code{Matching} package.}
}

Using the estimated propensity score, the individual observations in the dataset are matched to each other creating a new dataset in which the covariate distributions are balanced in respect to the two groups defined by \code{variable}. A simple Aalen-Johansen estimator is then used to calculate the confounder-adjusted survival curves.

Simulation results showed that this specific implementation of this method is the least efficient method contained in this R-Package. While it does produce unbiased estimates, the variation in these estimates is very high. We strongly suggest using one of the other methods implemented here.

In contrast to all other methods, bootstrapped standard error and confidence intervals are generally not unbiased here and are therefore not allowed. See Abadie & Imbens (2008) for a detailed explanation.

}
\value{
Returns a \code{data.frame} object containing the estimated adjusted survival probabilities for some points in time for each level of \code{variable}.
}
\references{
Angela Winnett and Peter Sasieni. "Adjusted Nelson-Aalen Estimates with Retrospective Matching". In: Journal of the American Statistical Association 97.457 (2002), pp. 245-256.

Peter C. Austin, Neal Thomas, and Donald B. Rubin. "Covariate-Adjusted Survival Analyses in Propensity-Score Matched Samples: Imputing Potential Time-To-Event Outcomes". In: Statistical Methods in Medical Research 29.3 (2020), pp. 728-751.

Alberto Abadie and Guido W. Imbens. "On the Failure of the Bootstrap for Matching Estimators". In: Econometrica 76.6 (2008), pp. 1537–1557.
}
\author{
The wrapper function was written by Robin Denz, the \code{Match} function (which this wrapper is build around) was written by other people. See \code{?Match} for more details.
}

\seealso{
\code{\link[Matching]{Match}}, \code{\link[survival]{survfit}}
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
