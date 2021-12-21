\name{surv_iptw_km}
\alias{surv_iptw_km}

\title{
Inverse Probability of Treatment Weighted Kaplan-Meier estimates
}
\description{

This page explains the details of estimating inverse probability of treatment weighted survival curves using a weighted version of the Kaplan-Meier estimator for single event time-to-event data. The function is called internally when using \code{method="iptw_km"} in the \code{\link[=adjustedsurv]{adjustedsurv}} function. It is highly recommended to use \code{\link[=adjustedsurv]{adjustedsurv}} instead of calling this function directly.

When using \code{method="iptw_km"} the \code{treatment_model} argument has to be specified in the \code{adjustedsurv} call.

}
\usage{
surv_iptw_km(data, variable, ev_time, event, conf_int,
             conf_level=0.95, times=NULL, treatment_model,
             weight_method="ps", stabilize=TRUE,
             trim=FALSE, ...)
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
A logical variable, indicating whether the asymptotic variances and confidence intervals of the survival probabilities should be calculated.
  }
  \item{conf_level}{
A number specifying the confidence level of asymptotic confidence intervals.
  }
  \item{times}{
A numeric vector of time points for which to return the survival probability estimates or \code{NULL} (default). If \code{NULL} the whole curve will be estimated according to the normal procedure. When custom time points are specified the whole curve is estimated anyways and the time specific estimates are simply read off that curve.
}
  \item{treatment_model}{
Must be either a model object with \code{variable} as response variable, a vector of weights or a formula which can be passed to \code{WeightIt}.
}
  \item{weight_method}{
Method used in \code{WeightIt} function call. Ignored if \code{treatment_model} is not a formula object.
}
  \item{stabilize}{
Whether to stabilize the weights or not. Is set to \code{TRUE} by default. Stabilizing weights ensures that the sum of all weights is equal to the original sample size. It has no effect on point estimates, only on the asymptotic variance calculations and confidence intervals. Stabilized weights are calculated using \code{w * length(w) / sum(w)} where \code{w} is a vector of weights.
  }
  \item{trim}{
Can be either \code{FALSE} (default) or a numeric value at which to trim the weights. If \code{FALSE}, weights are used as calculated or supplied. If a numeric value is supplied, all weights that are bigger than \code{trim} are set to \code{trim} before the analysis is carried out. Useful when some weights are extremely large.
  }
  \item{...}{
Further arguments passed to \code{\link[WeightIt]{weightit}}.
  }
}
\details{

\itemize{
\item{\strong{Type of Adjustment:} Requires a model describing the treatment assignment mechanism. This must be either a \code{\link{glm}} or \code{\link{multinom}} object.}
\item{\strong{Doubly-Robust:} Estimates are not Doubly-Robust.}
\item{\strong{Categorical groups:} Any number of levels in \code{variable} are allowed. Must be a factor variable.}
\item{\strong{Approximate Variance:} Calculations to approximate the variance and confidence intervals are available.}
\item{\strong{Allowed Time Values:} Allows both continuous and integer time.}
\item{\strong{Bounded Estimates:} Estimates are guaranteed to be bounded in the 0 to 1 probability range.}
\item{\strong{Monotone Function:} Estimates are guaranteed to be monotone.}
\item{\strong{Dependencies:} This method does not depend on other packages directly. However the \pkg{WeightIt} package is required if \code{treatment_model} is a formula object.}
}

This method works by modelling the treatment assignment mechanism. Adjusted survival curves are calculated by first estimating appropriate case-weights for each observation in \code{data}. This can be done using inverse probability of treatment weights using the propensity score (usually estimated using a logistic regression model) or by some other method (see \code{?weightit}). Those weights are used in a weighted version of the Kaplan-Meier estimator proposed by Xie and Liu (2005). If the weights are correctly estimated the resulting estimates will be unbiased. The only difference to the \code{\link[=surv_iptw_cox]{iptw_cox}} method is a slightly different weighting approach.

Asymptotic variances are calculated using the equations given in Xie and Liu (2005). It is also recommended to use stabilized weights by using \code{stabilize=TRUE} (the default value). More information can be found in the cited literature.
}
\value{
Returns a \code{list} object containing a \code{data.frame} with the estimated adjusted survival probabilities for some points in time for each level of \code{variable} and the weights used in the analysis.
}
\references{
Jun Xie and Chaofeng Liu (2005). "Adjusted Kaplan-Meier Estimator and Log- Rank Test with Inverse Probability of Treatment Weighting for Survival Data". In: Statistics in Medicine 24, pp. 3089-3110

Stanley Xu, Colleen Ross abd Marsha A. Raebel, Susan Shetterly, Christopher Blanchette, and David Smith (2010). "Use of Stabilized Inverse Propensity Scores as Weights to Directly Estimate Relative Risk and Its Confidence Intervals". In: Value in Health 13.2, pp. 273-277
}
\author{
Robin Denz
}

\seealso{
\code{\link[WeightIt]{weightit}}
}
\examples{
\dontrun{
# simulate some data as example
sim_dat <- sim_confounded_surv(n=500, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# estimate a treatment assignment model
glm_mod <- glm(group ~ x1 + x3 + x5 + x6, data=sim_dat, family="binomial")

# use it to calculate adjusted survival curves
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="iptw_km",
                        treatment_model=glm_mod)

# Alternatively, use custom weights
# In this example we use weights calculated using the propensity score,
# which is equal to using the glm model directly in the function
ps_score <- glm_mod$fitted.values
weights <- ifelse(sim_dat$group==1, 1/ps_score, 1/(1-ps_score))

adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="iptw_km",
                        treatment_model=weights)

# And a third alternative: use the WeightIt package
# here an example with equal results to the ones above:
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="iptw_km",
                        treatment_model=group ~ x1 + x3 + x5 + x6,
                        weight_method="ps")

# here an example using Optimization-Based Weighting:
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="iptw_km",
                        treatment_model=group ~ x1 + x3 + x5 + x6,
                        weight_method="optweight")
}
}
