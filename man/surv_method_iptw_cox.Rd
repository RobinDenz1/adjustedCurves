\name{surv_method_iptw_cox}
\alias{surv_method_iptw_cox}

\title{
Inverse Probability of Treatment Weighted Survival using Cox-Regression
}
\description{
This method works by modelling the treatment assignment mechanism. The individual observations are weighted using inverse probability of treatment weights or other balancing weights (see details) and a univariate cox-regression model. If the model is correctly specified, unbiased results will be obtained. The only difference to the \code{surv_method_iptw_km} is a slightly different weighting approach. They are essentially equal in large samples. See details for more information.
}
\usage{
surv_method_iptw_cox(data, variable, ev_time, event, conf_int,
                     conf_level=0.95, treatment_model, weight_method="ps",
                     stabilize=T, ...)
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
A logical variable, indicating whether the asymptotic variances and confidence intervals of the survival probabilities should be calculated. Not available for all methods.
  }
  \item{conf_level}{
A number specifying the confidence level of asymptotic and/or bootstrap confidence intervals.
  }
  \item{treatment_model}{
Must be either a model object with \code{variable} as response variable, a vector of weights or a formula which can be passed to \code{WeightIt}.
}
  \item{weight_method}{
Method used in \code{WeightIt} function call. Ignored if \code{treatment_model} is not a formula object.
}
  \item{stabilize}{
Whether to stabilize the weights or not. Is set to \code{TRUE} by default. Stabilizing has no effect on the estimates, only on the asymptotic variance calculations and confidence intervals. Without stabilizing, those estimates are biased downward.
  }
  \item{...}{
Further arguments passed to \code{WeightIt}.
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
\item{\strong{Dependencies:} This method relies on the \code{\link{survival}} package. Additionally, the \code{WeightIt} package is required if \code{treatment_model} is a formula object.}
}

Adjusted survival curves are calculated by first estimating appropriate case-weights for each observation in \code{data}. This can be done using inverse probability of treatment weights using the propensity score (usually estimated using a logistic regression model) or by some other method (see \code{?weightit}). Those estimates are then used to fit a weighted Cox-Regression model, using only \code{variable} as independent variable. Survival Curves based on this model are estimated using the method implemented in the \code{survfit.coxph} function. More information can be found in the literature listed under "references".

By default this method uses a a robust sandwich-type variance estimator (\code{robust=TRUE} in the \code{coxph} function call) to calculate the standard error used in the construction of confidence intervals. This estimator has been shown to be biased by Austin (2016). Coupled with stabilized weights however (\code{stabilize=TRUE}) this gives conservative estimates of the variance and confidence intervals (Xu et al. 2010). It is still recommended to use bootstrap confidence intervals instead. This can be done by setting \code{bootstrap=TRUE} in the \code{adjustedsurv} function call.
}
\value{
Returns a \code{data.frame} object containing the estimated adjusted survival probabilities for some points in time for each level of \code{variable}.
}
\references{
Stephen R. Cole and Miguel A. Hernan. "Adjusted Survival Curves with Inverse Probability Weights". In: Computer Methods and Programs in Biomedicine 2003.75 (2004), pp. 45-49.

Peter C. Austin. "Variance Estimation when Using Inverse Probabiltiy of Treatment Weighting (IPTW) with Survival Analysis". In: Statistics in Medicine 35 (2016), pp. 5642-5655.

Stanley Xu, Colleen Ross abd Marsha A. Raebel, Susan Shetterly, Christopher Blanchette, and David Smith. "Use of Stabilized Inverse Propensity Scores as Weights to Directly Estimate Relative Risk and Its Confidence Intervals". In: Value in Health 13.2 (2010), pp. 273-277.
}
\author{
Robin Denz
}

\seealso{
\code{\link{WeightIt}}, \code{\link{coxph}}, \code{\link{survfit.coxph}}
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
                        method="iptw_cox",
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
                        method="iptw_cox",
                        treatment_model=weights)

# And a third alternative: use the WeightIt package
# here an example with equal results to the ones above:
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="iptw_cox",
                        treatment_model=group ~ x1 + x3 + x5 + x6,
                        weight_method="ps")

# here an example using Optimization-Based Weighting:
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="iptw_cox",
                        treatment_model=group ~ x1 + x3 + x5 + x6,
                        weight_method="optweight")
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
