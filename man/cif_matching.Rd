\name{cif_matching}
\alias{cif_matching}

\title{
Using Matching to Calculate Adjusted CIFs
}
\description{

This page explains the details of estimating adjusted cumulative incidence functions using propensity-score matching in a competing risks setting. The function is called internally when using \code{method="matching"} in the \code{\link[=adjustedcif]{adjustedcif}} function. It is highly recommended to use \code{\link[=adjustedcif]{adjustedcif}} instead of calling this function directly.

}
\usage{
cif_matching(data, variable, ev_time, event, cause,
             conf_int, conf_level=0.95, times,
             treatment_model, gtol=0.001, ...)
}

\arguments{
  \item{data}{
A \code{data.frame} object containing the needed time-to-event data in standard format.
}
  \item{variable}{
A character string specifying the variable by which the CIFs should be grouped. Must be a valid column name of \code{data}. The variable in question can be either a dichotomous integer variable or a factor variable with exactly two levels.
}
  \item{ev_time}{
A character string specifying the variable indicating the time-to-event or time-to-censoring. Must be a valid column name of \code{data}.
}
  \item{event}{
A character string specifying the numeric event indicator. The censoring indicator should be coded as 0 and all other events of interest as 1, 2, etc. Must be a valid column name of \code{data}.
  }
  \item{cause}{
The cause of interest for which the cumulative incidence functions should be calculated. Should be a number that appears in the \code{event} column of \code{data}.
  }
  \item{conf_int}{
Currently not used because asymptotic variance calculations for this estimator have not been implemented.
  }
  \item{conf_level}{
A number specifying the confidence level of asymptotic confidence intervals.
  }
  \item{times}{
A numeric vector of time points for which to return the cumulative incidence estimates or \code{NULL} (default). If \code{NULL} the whole curve will be estimated according to the usual procedure. When custom time points are specified the whole curve is estimated anyways and the time specific estimates are simply read off that curve.
}
  \item{treatment_model}{
Must be either a model object with \code{variable} as response variable or a vector of previously estimated propensity scores.
}
  \item{gtol}{
Tolerance at which estimated treatment assignment probabilities are truncated. Every propensity score bigger than 1 - \code{gtol} is set to 1 - \code{gtol} and every propensity score smaller than \code{gtol} is set to \code{gtol}. Useful when there are extreme propensity scores close to 0 or 1.
  }
  \item{...}{
Further arguments passed to the \code{Match} function of the \pkg{Matching} Package.
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
\item{\strong{Dependencies:} This method relies on both the \pkg{Matching} and the \pkg{cmprsk} packages.}
}

Using the estimated propensity score, the individual observations in the dataset are matched to each other creating a new dataset in which the covariate distributions are balanced in respect to the two groups defined by \code{variable}. A simple Aalen-Johansen estimator is then used to calculate the cause-specific confounder-adjusted CIFs.

Simulation results showed that this specific implementation of this method is the least efficient method contained in this R-Package. While it does produce unbiased estimates, the variation in these estimates is very high. We strongly suggest using one of the other methods implemented here.

In contrast to all other methods, bootstrapped standard error and confidence intervals are generally not unbiased here and are therefore not allowed. See Abadie & Imbens (2008) for a detailed explanation.

}
\value{
Returns a \code{list} object containing a \code{data.frame} with estimated CIFs for some points in time for each level of \code{variable} and the object created using the \code{\link[Matching]{Match}} function.
}
\references{
Peter C. Austin and Jason P. Fine. "Propensity-Score Matching with Competing Risks in Survival Analysis". In: Statistics in Medicine 38 (2019), pp. 751–777.

Alberto Abadie and Guido W. Imbens. "On the Failure of the Bootstrap for Matching Estimators". In: Econometrica 76.6 (2008), pp. 1537–1557.
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
sim_dat <- sim_confounded_crisk(n=300, max_t=5)

# estimate treatment assignment model
glm_mod <- glm(group ~ x1 + x2 + x4 + x6, data=sim_dat, family="binomial")

# calculate adjusted CIFs
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      cause=1,
                      method="matching",
                      treatment_model=glm_mod)
plot(adjcif)

# Alternatively, supply the propensity score directly
# Here we use the logistic regression to calculate it, so we get
# exactly the same result. The propensity score can be calculated in
# any other way in practice, allowing flexibility
ps_score <- glm_mod$fitted.values

adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      cause=1,
                      method="matching",
                      treatment_model=ps_score)

# plot the curves
plot(adjcif)
}
}
