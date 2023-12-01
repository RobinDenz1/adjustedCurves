\name{adjusted_surv_quantile}
\alias{adjusted_surv_quantile}

\title{
Estimate Confounder-Adjusted Survival Time Quantiles
}
\description{
This function can be utilized to estimate confounder-adjusted survival time quantiles, including the median survival time, given previously estimated adjusted survival curves.
}
\usage{
adjusted_surv_quantile(adjsurv, p=0.5, conf_int=FALSE,
                       use_boot=FALSE, interpolation="steps")
}

\arguments{
  \item{adjsurv}{
An \code{adjustedsurv} object created using the \code{adjustedsurv} function.
  }
  \item{p}{
The quantile of interest. To calculate the median survival time, set this parameter to 0.5 (default). Multiple values in form of a numeric vector are allowed.
  }
  \item{conf_int}{
Whether to calculate confidence intervals or not. Those are calculated in the same way as the quantiles but using the confidence limit curves. This requires either that \code{conf_int=TRUE} or \code{bootstrap=TRUE} was used in the original \code{adjustedsurv} function call. Since this directly uses the previously estimated intervals, the same confidence level used in the original \code{adjustedsurv} call is used here.
  }
  \item{use_boot}{
Whether to use the bootstrap confidence interval estimates of the survival curves to estimate the confidence intervals of the survival time quantiles or not. Can only be used when \code{bootstrap=TRUE} was used in the original \code{adjustedsurv} function call. Ignored if \code{conf_int=FALSE}.
  }
  \item{interpolation}{
Either \code{"steps"} (default) or \code{"linear"}. This parameter controls how interpolation is performed. If this argument is set to \code{"steps"}, the curves will be treated as step functions. If it is set to \code{"linear"}, the curves wil be treated as if there are straight lines between the point estimates instead. Points that lie between estimated points will be interpolated accordingly.
  }
}
\details{

The median survival time is simply the time when half the patients are expected to be alive. The chance of surviving beyond that time is 50 percent. In general, any quantile can be calculated this way. Those can be read directly from the respective survival curve by drawing a straight line from the desired quantile \code{p} on the Y-Axis and reading the X-Axis value where this line intersects with the survival curve. The adjusted survival time quantile for group \eqn{z} (corresponding to the \code{variable} argument in the \code{\link{adjustedsurv}} function) is formally defined as:

\deqn{\hat{Q}_z(p) = min\left(t | \hat{S}_z(t) \leq p\right)}

where \eqn{\hat{S}_z(t)} is the estimated counterfactual survival function for \eqn{z}.

If the survival probability never drops below \code{p}, the survival time quantile cannot be calculated. This also applies to the confidence interval estimation. This function calculates this quantity automatically. When multiple imputation was used in the original function call, the survival time quantiles are read off the final pooled survival curves directly.

}
\value{

Returns a \code{data.frame} containing the columns \code{p} (the quantiles from the original function call), \code{group} (groups in \code{variable}) and \code{q_surv} (the survival time quantile).

If \code{conf_int=TRUE} was used it also includes the confidence limits in the \code{ci_lower} and \code{ci_upper} columns.

}
\references{
Omer Ben-Aharon, Racheli Magnezi, Moshe Leshno, and Daniel A. Goldstein (2019). "Median Survival or Mean Survival: Which Measure is the Most Appropriate for Patients, Physicians, and Policymakers?" In: The Oncologist 24, pp. 1469-1478

Zhongxue Chen and Guoyi Zhang (2016). "Comparing Survival Curves based on Medians". In: BMC Medical Research Methodology 16.33
}
\author{
Robin Denz
}

\seealso{
\code{\link{adjustedsurv}}
}
\examples{
library(adjustedCurves)
library(survival)

# simulate some data as example
set.seed(42)
sim_dat <- sim_confounded_surv(n=50, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# estimate a cox-regression for the outcome
cox_mod <- coxph(Surv(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6 + group,
                 data=sim_dat, x=TRUE)

# use it to calculate adjusted survival curves with bootstrapping
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="direct",
                        outcome_model=cox_mod,
                        conf_int=TRUE,
                        bootstrap=FALSE)

# calculate adjusted median survival times
adjusted_surv_quantile(adjsurv)

# calculate other quantiles + confidence intervals
adjusted_surv_quantile(adjsurv, conf_int=TRUE, p=c(0.2, 0.4))
}
