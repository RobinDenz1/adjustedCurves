\name{adjusted_rmst}
\alias{adjusted_rmst}

\title{
Calculate Confounder-Adjusted Restricted Mean Survival Times
}
\description{
This function can be utilized to calculate the confounder-adjusted restricted mean survival time, given previously estimated adjusted survival curves.
}
\usage{
adjusted_rmst(adjsurv, to, from=0, use_boot=F, conf_level=0.95)
}

\arguments{
  \item{adjsurv}{
An \code{adjustedsurv} object created using the \code{adjustedsurv} function.
  }
  \item{from}{
A number specifying the left side of the time interval of interest. See details.
  }
  \item{to}{
A number specifying the right side of the time interval of interest. See details.
  }
  \item{use_boot}{
Whether bootstrap estimates should be used to calculate the standard deviation and confidence intervals of the RMST estimates. Can only be used if \code{bootstrap=TRUE} was used in the \code{\link{adjustedsurv}} call.
  }
  \item{conf_level}{
A number specifying the confidence level of the bootstrap confidence intervals.
  }
}
\details{

The adjusted restricted mean survival times (RMST) are calculated by integrating the estimated adjusted survival curves in a specified interval. Let \eqn{Z} be the grouping variable (corresponding to the \code{variable} argument in the \code{\link{adjustedsurv}} function) with possible levels \eqn{Z = z \in \{0, 1, 2, ..., k\}}. \eqn{T} is defined as the time and \eqn{X} is a vector or matrix of covariates. \eqn{\hat{S}(T = t | Z = z, X)} denotes the estimated counterfactual survival function. The RMST is then defined as:

\deqn{RMST_{z} = \int_{from}^{to} \hat{S}(T = t | Z = z, X)}

It can be interpreted as the mean survival time of individuals in group \eqn{Z = z}. Note however that simply substracting the estimates from each other does not give a correct estimate of the area between the survival curves if the respective curves cross at some point. The \code{\link{test_curve_equality}} function can be used to calculate the actual area between the curves instead. See \code{?test_curve_equality} for more information.

Since the estimated survival functions are step functions the integrals are calculated exactly instead of relying on trapezoid approximations.

If the \code{adjsurv} object was created with \code{bootstrap = TRUE} in the \code{\link{adjustedsurv}} function, bootstrap confidence intervals and standard deviations for the RMSTs can be approximated by setting \code{use_boot} to \code{TRUE}. If bootstrap samples occur where the survival function is not estimated up to \code{to}, the bootstrap sample is discarded and not used in further calculations. Approximate variance calculations not relying on the bootstrap estimates are currently not implemented.

}
\value{
Returns a \code{adjusted_rmst} object containing the following object:

\item{rmsts}{A named vector containing the adjusted restricted mean survival time estimates.
}
\item{from}{\code{from} argument in function call.
}
\item{to}{\code{to} argument in function call.
}

If \code{use_boot} is set to \code{TRUE} and \code{bootstrap} was also set to \code{TRUE} in the original \code{adjustedsurv} function call, the following objects are also added to the output:

\item{conf_level}{The confidence level used when calculating the confidence intervals.
}
\item{n_boot}{The numbers of bootstrap iterations actually used to calculate the bootstrap estimates.
}
\item{booted_rmsts}{A \code{data.frame} containing the adjusted restricted mean survival time estimates for each bootstrap sample.
}
\item{rmsts_sd}{A named vector containing the bootstrapped standard deviations of the restricted mean survival times.
}
\item{rmsts_ci_lower}{Lower value of the bootstrap confidence intervals.
}
\item{rmsts_ci_upper}{Upper value of the bootstrap confidence intervals.
}

}
\references{
Sarah C. Conner, Lisa M. Sullivan, Emelia J. Benjamin, Michael P. LaValley, Sandro Galea, and Ludovic Trinquart. "Adjusted Restricted Mean Survival Times in Observational Studies". In: Statistics in Medicine 38 (2019), pp. 3832-3860.

Patrick Royston and Mahesh K. B. Parmar. "Restricted Mean Survival Time: An Alternative to the Hazard Ratio for the Design and Analysis of Randomized Trials with a Time-To-Event Outcome". In: BMC Medical Research Methodology 13.152 (2013).
}
\author{
Robin Denz
}

\seealso{
\code{\link{print.adjusted_rmst}}, \code{\link{adjustedsurv}}
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
adjrmst <- adjusted_rmst(adjsurv, from=0, to=1, use_boot=FALSE)
print(adjrmst)

# calculate adjusted restricted mean survival times from 0 to 1,
# including standard deviations and confidence intervals
adjrmst <- adjusted_rmst(adjsurv, from=0, to=1, use_boot=TRUE,
                         conf_level=0.95)
print(adjrmst)
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
