\name{cif_matching}
\alias{cif_matching}

\title{
Using Propensity-Score Matching to Calculate Adjusted CIFs
}
\description{

This page explains the details of estimating adjusted cumulative incidence functions using propensity-score matching in a competing risks setting (\code{method="matching"} in the \code{\link[=adjustedcif]{adjustedcif}} function). All regular arguments of the \code{adjustedcif} function can be used. Additionally, the \code{treatment_model} argument has to be specified in the \code{adjustedcif} call. Further arguments specific to this method are listed below.

}

\arguments{
  \item{treatment_model}{
[\strong{required}] Must be either a model object with \code{variable} as response variable or a vector of previously estimated propensity scores.
}
  \item{gtol}{
Tolerance at which estimated treatment assignment probabilities are truncated. Every propensity score bigger than 1 - \code{gtol} is set to 1 - \code{gtol} and every propensity score smaller than \code{gtol} is set to \code{gtol}. Useful when there are extreme propensity scores close to 0 or 1. Defaults to 0.001,
  }
  \item{...}{
Further arguments passed to the \code{Match} function of the \pkg{Matching} Package.
  }
}
\details{

\itemize{
\item{\strong{Type of Adjustment:} Requires a model describing the treatment assignment mechanism. This must be either a \code{\link{glm}} object or a vector of propensity scores.}
\item{\strong{Doubly-Robust:} Estimates are not Doubly-Robust.}
\item{\strong{Categorical groups:} Only two groups in \code{variable} are allowed. Must be a factor variable with exactly two levels.}
\item{\strong{Approximate Variance:} Calculations to approximate the variance and confidence intervals are currently not available. Bootstrapping can still be used to estimate the confidence intervals (see \code{?adjustedcif}).}
\item{\strong{Allowed Time Values:} Allows both continuous and integer time.}
\item{\strong{Bounded Estimates:} Estimates are guaranteed to be bounded in the 0 to 1 probability range.}
\item{\strong{Monotone Function:} Estimates are guaranteed to be monotone.}
\item{\strong{Dependencies:} This method relies on both the \pkg{Matching} and the \pkg{cmprsk} packages.}
}

Using the estimated propensity score, the individual observations in the dataset are matched to each other creating a new dataset in which the covariate distributions are balanced in respect to the two groups defined by \code{variable}. A simple Aalen-Johansen estimator is then used to calculate the confounder-adjusted CIFs. This corresponds to the method described in Austin & Fine (2019). Details on the algorithm used for matching can be found in the documentation of the \pkg{Matching} package.

Simulation results showed that this specific implementation of this method is the least efficient method contained in this R-Package. While it does produce unbiased estimates, the variation in these estimates is very high. We strongly suggest using one of the other methods implemented here.

}
\value{
Adds the following additional objects to the output of the \code{adjustedcif} function:

\itemize{
\item \code{match_object}: The object creates using the \code{Match} function.
\item \code{cuminc_object}: The \code{cuminc} object fit on the matched data.
}
}
\references{
Peter C. Austin and Jason P. Fine (2019). "Propensity-Score Matching with Competing Risks in Survival Analysis". In: Statistics in Medicine 38, pp. 751-777
}
\author{
Robin Denz
}

\seealso{
\code{\link[Matching]{Match}}, \code{\link[cmprsk]{cuminc}}
}
\examples{
library(adjustedCurves)
library(survival)

if (requireNamespace("Matching")) {

library(Matching)

set.seed(42)

# simulate some data as example
sim_dat <- sim_confounded_crisk(n=50, max_t=5)
sim_dat$group <- as.factor(sim_dat$group)

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
