\name{adjusted_rmst}
\alias{adjusted_rmst}

\title{
Calculate Confounder-Adjusted Restricted Mean Survival Times
}
\description{
This function can be utilized to calculate the confounder-adjusted restricted mean survival time, given previously estimated adjusted survival curves.
}
\usage{
adjusted_rmst(adjsurv, to, from=0, conf_int=FALSE,
              conf_level=0.95, interpolation="steps",
              subdivisions=1000)
}

\arguments{
  \item{adjsurv}{
An \code{adjustedsurv} object created using the \code{adjustedsurv} function.
  }
  \item{from}{
A number specifying the left side of the time interval of interest. See details. Usually this should be kept at 0 (default) to calculate the standard RMST. Should only be changed if there are good reasons for it.
  }
  \item{to}{
A number specifying the right side of the time interval of interest. See details.
  }
  \item{conf_int}{
Whether bootstrap estimates should be used to calculate the standard errors and confidence intervals of the RMST estimates. Can only be used if \code{bootstrap=TRUE} was used in the \code{\link{adjustedsurv}} call.
  }
  \item{conf_level}{
A number specifying the confidence level of the bootstrap confidence intervals.
  }
  \item{interpolation}{
Either \code{"steps"} (default) or \code{"linear"}. This parameter controls how interpolation is performed. If this argument is set to \code{"steps"}, the curves will be treated as step functions. If it is set to \code{"linear"}, the curves wil be treated as if there are straight lines between the point estimates instead. Points that lie between estimated points will be interpolated accordingly. Should usually be kept at \code{"steps"}. See Details.
  }
  \item{subdivisions}{
When using \code{interpolation="linear"}, the composite trapezoid method is used to calculate the integral of the curves. This argument defines the number of subdivisions that should be used when calculating the integral. Bigger number give a more accurate estimate, but take longer to compute.
  }
}
\details{

The adjusted restricted mean survival times (RMST) are calculated by integrating the estimated adjusted survival curves in a specified interval. Let \eqn{Z} be the grouping variable (corresponding to the \code{variable} argument in the \code{\link{adjustedsurv}} function) with possible levels \eqn{Z \in \{0, 1, 2, ..., k\}}. \eqn{T} is defined as the time and \eqn{\hat{S}_z(t)} denotes the estimated counterfactual survival function. The RMST is then defined as:

\deqn{RMST_{z}(to) = \int_{from=0}^{to} \hat{S}_z(t)dt}

It can be interpreted as the mean survival time of individuals in group \eqn{Z = z} in the interval [\code{from}, \code{to}]. Note however that simply subtracting the estimates from each other does not give a correct estimate of the area between the survival curves if the respective curves cross at some point. The \code{\link{adjusted_curve_test}} function can be used to calculate the actual area between the curves instead. See \code{?adjusted_curve_test} for more information.

\strong{\emph{Confidence Intervals}}

If the \code{adjsurv} object was created with \code{bootstrap=TRUE} in the \code{\link{adjustedsurv}} function, bootstrap confidence intervals and standard errors for the RMSTs can be approximated by setting \code{conf_int} to \code{TRUE}. If bootstrap samples occur where the survival function is not estimated up to \code{to}, the bootstrap sample is discarded and not used in further calculations. Approximate variance calculations not relying on the bootstrap estimates are currently not implemented.

\strong{\emph{Multiple Imputation}}

If multiple imputation was used when creating the \code{adjsurv} object, the analysis is carried out on all multiply imputed datasets and pooled using Rubins Rule. When bootstrapping was carried out as well, the pooled standard error over all imputed datasets is used in combination with the normal approximation to re-calculate the bootstrap confidence intervals.

\strong{\emph{Competing Risks}}

This function cannot be used with \code{adjustedcif} objects, because the survival probability cannot be estimated in an unbiased way when competing risks are present. However, a very similar quantity, the \emph{adjusted restricted mean time lost}, can be calculated using the \code{\link{adjusted_rmtl}} function.

\strong{\emph{Graphical Displays}}

A plot of the RMST over time (with changing values for the \code{to} argument) can be produced using the \code{\link{plot_rmst_curve}} function.

\strong{\emph{Computational Details}}

Since the estimated survival functions are usually step functions, the integrals are by default calculated exactly instead of relying on trapezoid approximations. By setting \code{interpolation="linear"} the trapezoid version can be calculated instead. This, however, should only be used when the curves are not step functions. For example, this is the case when using a parametric survival model with \code{method="direct"}. The produced curves are smooth functions then, requiring linear interpolation for correct integral calculations.

}
\value{
Returns a \code{data.frame} containing the columns \code{group} (groups in \code{variable}) and \code{rmst} (the estimated restricted mean survival time).

If \code{conf_int=TRUE} was used it additionally contains the columns \code{se} (the standard error of the restricted mean survival time), \code{ci_lower} (lower limit of the confidence interval), \code{ci_upper} (upper limit of the confidence interval) and \code{n_boot} (the actual number of bootstrap estimates used).
}
\references{
Sarah C. Conner, Lisa M. Sullivan, Emelia J. Benjamin, Michael P. LaValley, Sandro Galea, and Ludovic Trinquart (2019). "Adjusted Restricted Mean Survival Times in Observational Studies". In: Statistics in Medicine 38, pp. 3832-3860

Patrick Royston and Mahesh K. B. Parmar (2013). "Restricted Mean Survival Time: An Alternative to the Hazard Ratio for the Design and Analysis of Randomized Trials with a Time-To-Event Outcome". In: BMC Medical Research Methodology 13.152
}
\author{
Robin Denz
}

\seealso{
\code{\link{adjustedsurv}}, \code{\link{plot_rmst_curve}}
}
\examples{
\dontrun{
library(survival)

# simulate some data as example
sim_dat <- sim_confounded_surv(n=500, max_t=1.2)
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
                        n_boot=500)

# calculate adjusted restricted mean survival times from 0 to 1
adjrmst <- adjusted_rmst(adjsurv, from=0, to=1, conf_int=FALSE)

# calculate adjusted restricted mean survival times from 0 to 1,
# including standard errors and confidence intervals
adjrmst <- adjusted_rmst(adjsurv, from=0, to=1, conf_int=TRUE,
                         conf_level=0.95)
}
}
