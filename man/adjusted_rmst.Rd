\name{adjusted_rmst}
\alias{adjusted_rmst}

\title{
Estimate Confounder-Adjusted Restricted Mean Survival Times
}
\description{
This function can be utilized to estimate the confounder-adjusted restricted mean survival time, given previously estimated adjusted survival curves.
}
\usage{
adjusted_rmst(adjsurv, to, from=0, conf_int=FALSE,
              conf_level=0.95, interpolation="steps",
              contrast="none", group_1=NULL,
              group_2=NULL)
}
\arguments{
  \item{adjsurv}{
An \code{adjustedsurv} object created using the \code{adjustedsurv} function.
  }
  \item{from}{
A single number specifying the left side of the time interval of interest. See details. Usually this should be kept at 0 (default) to estimate the standard RMST. Should only be changed if there are good reasons for it.
  }
  \item{to}{
One or multiple numbers specifying the right side of the time interval of interest. If a vector of numbers is supplied, the adjusted RMST will be estimated for each value of \code{to}. See details.
  }
  \item{conf_int}{
Whether bootstrap estimates should be used to estimate the standard errors and confidence intervals of the RMST estimates. Can only be used if \code{bootstrap=TRUE} was used in the \code{\link{adjustedsurv}} call.
  }
  \item{conf_level}{
A number specifying the confidence level of the bootstrap confidence intervals.
  }
  \item{interpolation}{
Either \code{"steps"} (default) or \code{"linear"}. This parameter controls how interpolation is performed. If this argument is set to \code{"steps"}, the curves will be treated as step functions. If it is set to \code{"linear"}, the curves wil be treated as if there are straight lines between the point estimates instead. Points that lie between estimated points will be interpolated accordingly. Should usually be kept at \code{"steps"}. See Details.
  }
  \item{contrast}{
A single character string, specifying which contrast should be estimated. Needs to be one of \code{"none"} (estimate no contrasts, just return the adjusted RMST, the default), \code{"diff"} (estimate the difference) or \code{"ratio"} (estimate the ratio). When \code{conf_int=TRUE} is also specified and bootstrapping was performed in the original \code{adjustedsurv} call, this function will also estimate the corresponding standard error, the confidence interval and a p-value testing whether the difference is equal to 0 (or the ratio is equal to 1). To specify which difference/ratio should be calculated, the \code{group_1} and \code{group_2} arguments can be used. By default, the difference/ratio between the first and second level in \code{variable} is computed.
  }
  \item{group_1}{
Optional argument to get a specific difference or ratio. This argument takes a single character string specifying one of the levels of the \code{variable} used in the original \code{adjustedsurv} or \code{adjustedcif} function call. This group will be subtracted from. For example if \code{group_1="A"} and \code{group_2="B"} and \code{contrast="diff"} the difference \code{A - B} will be used. If \code{NULL}, the order of the factor levels in the original \code{data} determines the order. Ignored if \code{contrast="none"}.
  }
  \item{group_2}{
Also a single character string specifying one of the levels of \code{variable}. This corresponds to the right side of the difference/ratio equation. See argument \code{group_2}. Ignored if \code{contrast="none"}.
  }
}
\details{
The adjusted restricted mean survival times (RMST) are estimated by integrating the estimated adjusted survival curves in a specified interval. Let \eqn{Z} be the grouping variable (corresponding to the \code{variable} argument in the \code{\link{adjustedsurv}} function) with possible levels \eqn{Z \in \{0, 1, 2, ..., k\}}. \eqn{T} is defined as the time and \eqn{\hat{S}_z(t)} denotes the estimated counterfactual survival function. The RMST is then defined as:

\deqn{RMST_{z}(to) = \int_{from=0}^{to} \hat{S}_z(t)dt}

It can be interpreted as the mean survival time of individuals in group \eqn{Z = z} in the interval [\code{from}, \code{to}]. Note however that simply subtracting the estimates from each other does not give a correct estimate of the area between the survival curves if the respective curves cross at some point. The \code{\link{adjusted_curve_test}} function can be used to calculate the actual area between the curves instead. See \code{?adjusted_curve_test} for more information.

\strong{\emph{Confidence Intervals}}

If the \code{adjsurv} object was created with \code{bootstrap=TRUE} in the \code{\link{adjustedsurv}} function, bootstrap confidence intervals and standard errors for the RMSTs can be approximated by setting \code{conf_int} to \code{TRUE}. If bootstrap samples occur where the survival function is not estimated up to \code{to}, the bootstrap sample is discarded and not used in further calculations. Approximate variance calculations not relying on the bootstrap estimates are currently not implemented. When using \code{contrast="diff"} the standard error of the difference between the two RMST values is approximated by \eqn{SE_{group_1 - group_2} = \sqrt{SE_{group_1}^2 + SE_{group_2}^2}}. When using \code{contrast="ratio"} the confidence intervals are calculated using the approximate formula given by Fieller (1954), assuming that the values are independent.

\strong{\emph{More than Two Groups}}

If more than two groups are present in \code{variable}, all other comparisons except for \code{group_1 vs. group_2} are ignored. If multiple comparisons are desired, the user needs to call this function multiple times and adjust the \code{group_1} and \code{group_2} arguments accordingly.

\strong{\emph{Multiple Imputation}}

If multiple imputation was used when creating the \code{adjsurv} object, the analysis is carried out on all multiply imputed datasets and pooled using Rubins Rule. When bootstrapping was carried out as well, the pooled standard error over all imputed datasets is used in combination with the normal approximation to re-calculate the bootstrap confidence intervals.

\strong{\emph{Competing Risks}}

This function cannot be used with \code{adjustedcif} objects, because the survival probability cannot be estimated in an unbiased way when competing risks are present. However, a very similar quantity, the \emph{adjusted restricted mean time lost}, can be calculated using the \code{\link{adjusted_rmtl}} function.

\strong{\emph{Graphical Displays}}

A plot of the RMST over time (with changing values for the \code{to} argument) can be produced using the \code{\link{plot_rmst_curve}} function.

\strong{\emph{Computational Details}}

Instead of relying on numerical integration, this function uses exact calculations. This is achieved by using either step-function interpolation (\code{interpolation="steps"}, the default) or linear interpolation (\code{interpolation="linear"}). In the former case, the integral is simply the sum of the area of the squares defined by the step function. In the second case, the integral is simply the sum of the area of the rectangles. Either way, there is no need for approximations. In some situations (for example when using parametric survival models with \code{method="direct"}), the curves are not step functions. In this case the \code{interpolation} argument should be set to \code{"linear"}.

}
\value{
Returns a \code{data.frame} containing the columns \code{to} (the supplied \code{to} values), \code{group} (groups in \code{variable}) and \code{rmst} (the estimated restricted mean survival time).

If \code{conf_int=TRUE} was used it additionally contains the columns \code{se} (the standard error of the restricted mean survival time), \code{ci_lower} (lower limit of the confidence interval), \code{ci_upper} (upper limit of the confidence interval) and \code{n_boot} (the actual number of bootstrap estimates used).

If \code{contrast="diff"} was used, it instead returns a \code{data.frame} that contains the columns \code{to}, \code{diff} (the difference between the RMST values), \code{se} (the standard error of the difference), \code{ci_lower} (lower limit of the confidence interval), \code{ci_upper} (upper limit of the confidence interval) and \code{p_value} (the p-value of the one-sample t-test). The same results are presented when using \code{contrast="ratio"}, except that the \code{diff} column is named \code{ratio} and that there is no \code{se} column.
}
\references{
Sarah C. Conner, Lisa M. Sullivan, Emelia J. Benjamin, Michael P. LaValley, Sandro Galea, and Ludovic Trinquart (2019). "Adjusted Restricted Mean Survival Times in Observational Studies". In: Statistics in Medicine 38, pp. 3832-3860

Patrick Royston and Mahesh K. B. Parmar (2013). "Restricted Mean Survival Time: An Alternative to the Hazard Ratio for the Design and Analysis of Randomized Trials with a Time-To-Event Outcome". In: BMC Medical Research Methodology 13.152

Edgar C. Fieller (1954). "Some Problems in Interval Estimation". In: Journal of the Royal Statistical Society, Series B 16.2, pp. 175-185
}
\author{
Robin Denz
}

\seealso{
\code{\link{adjustedsurv}}, \code{\link{plot_rmst_curve}}
}
\examples{
library(adjustedCurves)
library(survival)

set.seed(42)

# simulate some data as example
sim_dat <- sim_confounded_surv(n=30, max_t=1.2)
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
                        conf_int=FALSE,
                        bootstrap=TRUE,
                        n_boot=10) # n_boot should be much higher in reality

# calculate adjusted restricted mean survival times from 0 to 0.5
adjrmst <- adjusted_rmst(adjsurv, from=0, to=0.5, conf_int=FALSE)

# calculate adjusted restricted mean survival times from 0 to 0.5
# and from 0 to 0.7 simultaneously
adjrmst <- adjusted_rmst(adjsurv, from=0, to=c(0.5, 0.7), conf_int=FALSE)

# calculate adjusted restricted mean survival times from 0 to 0.5,
# including standard errors and confidence intervals
adjrmst <- adjusted_rmst(adjsurv, from=0, to=0.5, conf_int=TRUE,
                         conf_level=0.95)

# calculate difference between adjusted restricted mean survival times
# from 0 to 0.5 in the two groups
adjrmst <- adjusted_rmst(adjsurv, from=0, to=0.5, conf_int=FALSE,
                         contrast="diff")

# calculate ratio between adjusted restricted mean survival times
# from 0 to 0.5 in the two groups
adjrmst <- adjusted_rmst(adjsurv, from=0, to=0.5, conf_int=FALSE,
                         contrast="ratio")
}
