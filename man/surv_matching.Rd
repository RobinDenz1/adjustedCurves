\name{surv_matching}
\alias{surv_matching}

\title{
Using Propensity-Score Matching to Calculate Adjusted Survival Curves
}
\description{

This page explains the details of estimating adjusted survival curves using propensity-score matching for single event time-to-event data (\code{method="matching"} in the \code{\link[=adjustedsurv]{adjustedsurv}} function). All regular arguments of the \code{adjustedsurv} function can be used. Additionally, the \code{treatment_model} argument has to be specified in the \code{adjustedsurv} call. Further arguments specific to this method are listed below.

}

\arguments{
  \item{treatment_model}{
[\strong{required}] Must be either a model object with \code{variable} as response variable or a vector of previously estimated propensity scores.
}
  \item{gtol}{
Tolerance at which estimated treatment assignment probabilities are truncated. Every propensity score bigger than 1 - \code{gtol} is set to 1 - \code{gtol} and every propensity score smaller than \code{gtol} is set to \code{gtol}. Useful when there are extreme propensity scores close to 0 or 1. Defaults to 0.001.
  }
  \item{...}{
Further arguments passed to the \code{\link[Matching]{Match}} function of the \pkg{Matching} Package.
  }
}
\details{

\itemize{
\item{\strong{Type of Adjustment:} Requires a model describing the treatment assignment mechanism. This must be either a \code{\link{glm}} object.}
\item{\strong{Doubly-Robust:} Estimates are not Doubly-Robust.}
\item{\strong{Categorical groups:} Only two groups in \code{variable} are allowed. Must be a factor variable with exactly two levels.}
\item{\strong{Approximate Variance:} Calculations to approximate the variance and confidence intervals are currently not available. Bootstrap confidence intervals can however be calculated with all supported models. See \code{?adjustedsurv} for more information on bootstrapping.}
\item{\strong{Allowed Time Values:} Allows both continuous and integer time.}
\item{\strong{Bounded Estimates:} Estimates are guaranteed to be bounded in the 0 to 1 probability range.}
\item{\strong{Monotone Function:} Estimates are guaranteed to be monotone.}
\item{\strong{Dependencies:} This method relies on the \pkg{Matching} package.}
}

Using the estimated propensity score, the individual observations in the dataset are matched to each other creating a new dataset in which the covariate distributions are balanced in respect to the two groups defined by \code{variable}. A simple Kaplan-Meier estimator is then used to calculate the confounder-adjusted survival curves. This corresponds to the method described in Austin (2014). Details on the algorithm used for matching can be found in the documentation of the \pkg{Matching} package.

We choose not to implement other matching based estimators (see Winnett & Sasieni (2002), Galimberti et al. (2002) and Austin (2020)) because of the wide range of matching algorithms and parameters. Trying to automate the matching process in a function like this would, in our opinion, disrupt the workflow of the user while also encouraging suboptimal practices. We however included this simple version of a matching estimator as a reference and to raise the awareness that using matching is a valid method to obtain adjusted survival curves.

Simulation studies have shown that this particular method as implemented here is significantly less efficient than other methods included in this R-Package. While it does produce unbiased estimates, the variation in these estimates is very high. We suggest using one of the other available methods.

}
\value{
Adds the following additional objects to the output of the \code{adjustedsurv} function:

\itemize{
\item \code{match_object}: The object creates using the \code{Match} function.
\item \code{survfit_object}: The \code{survfit} object fit on the matched data.
}
}
\references{
Peter C. Austin (2014). "The Use of Propensity Score Methods with Survival or Time-To-Event Outcomes: Reporting Measures of Effect Similar to those Used in Randomized Experiments". In: Statistics in Medicine 33, pp. 1242-1258

Angela Winnett and Peter Sasieni (2002). "Adjusted Nelson-Aalen Estimates with Retrospective Matching". In: Journal of the American Statistical Association 97.457, pp. 245-256

Stefania Galimberti, Peter Sasieni, and Maria Grazia Valsecchi (2002). "A Weighted Kaplan-Meier Estimator for Matched Data with Application to the Comparison of Chemotherapy and Bone-Marrow Transplant in Leukaemia". In: Statistics in Medicine 21, pp. 3847-3864

Peter C. Austin, Neal Thomas, and Donald B. Rubin (2020). "Covariate-Adjusted Survival Analyses in Propensity-Score Matched Samples: Imputing Potential Time- To-Event Outcomes". In: Statistical Methods in Medical Research 29.3, pp. 728-751
}
\author{
Robin Denz
}

\seealso{
\code{\link[Matching]{Match}}, \code{\link[survival]{survfit}}
}
\examples{
library(adjustedCurves)

if (requireNamespace("Matching")) {

library(Matching)

set.seed(42)

# simulate some data as example
sim_dat <- sim_confounded_surv(n=50, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

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
