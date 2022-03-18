\name{adjusted_rmst}
\alias{adjusted_rmst}

\title{
Calculate Confounder-Adjusted Restricted Mean Survival Times
}
\description{
This function can be utilized to calculate the confounder-adjusted restricted mean survival time, given previously estimated adjusted survival curves.
}
\usage{
adjusted_rmst(adjsurv, to, from=0, use_boot=FALSE,
              conf_level=0.95)
}

\arguments{
  \item{adjsurv}{
An \code{adjustedsurv} object created using the \code{adjustedsurv} function.
  }
  \item{from}{
A number specifying the left side of the time interval of interest. See details. Usually this should be kept at 0 (default) unless there are good reasons to change it.
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

The adjusted restricted mean survival times (RMST) are calculated by integrating the estimated adjusted survival curves in a specified interval. Let \eqn{Z} be the grouping variable (corresponding to the \code{variable} argument in the \code{\link{adjustedsurv}} function) with possible levels \eqn{Z \in \{0, 1, 2, ..., k\}}. \eqn{T} is defined as the time and \eqn{\hat{S}_z(t)} denotes the estimated counterfactual survival function. The RMST is then defined as:

\deqn{RMST_{z} = \int_{from}^{to} \hat{S}_z(t)}

It can be interpreted as the mean survival time of individuals in group \eqn{Z = z} in the interval [\code{from}, \code{to}]. Note however that simply subtracting the estimates from each other does not give a correct estimate of the area between the survival curves if the respective curves cross at some point. The \code{\link{adjusted_curve_diff}} function can be used to calculate the actual area between the curves instead. See \code{?adjusted_curve_diff} for more information.

The associated \code{print} and \code{summary} methods can be used to obtain a neat data.frame of the most important quantities. There is also a small \code{plot} method to graphically display the results.

\strong{\emph{Confidence Intervals}}

If the \code{adjsurv} object was created with \code{bootstrap = TRUE} in the \code{\link{adjustedsurv}} function, bootstrap confidence intervals and standard errors for the RMSTs can be approximated by setting \code{use_boot} to \code{TRUE}. If bootstrap samples occur where the survival function is not estimated up to \code{to}, the bootstrap sample is discarded and not used in further calculations. Approximate variance calculations not relying on the bootstrap estimates are currently not implemented.

\strong{\emph{Multiple Imputation}}

If multiple imputation was used when creating the \code{adjsurv} object, the analysis is carries out on all multiply imputed datasets and pooled using Rubins Rule. When bootstrapping was carried out as well, the pooled standard error over all imputed datasets is used in combination with the normal approximation to re-calculate the bootstrap confidence intervals.

\strong{\emph{Competing Risks}}

This function cannot be used with \code{adjustedcif} objects, because the survival probability cannot be estimated in an unbiased way when competing risks are present. However, a very similar quantity, the \emph{adjusted restricted mean time lost}, can be calculated using the \code{\link{adjusted_rmtl}} function.

\strong{\emph{Computational Details}}

Since the estimated survival functions are step functions the integrals are calculated exactly instead of relying on trapezoid approximations.

}
\value{
Returns a \code{adjusted_rmst} object containing the following object:

\item{auc}{A named vector containing the adjusted restricted mean survival time estimates.
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
\item{booted_auc}{A \code{data.frame} containing the adjusted restricted mean survival time estimates for each bootstrap sample.
}
\item{auc_se}{A named vector containing the bootstrapped standard errors of the restricted mean survival times.
}
\item{auc_ci_lower}{Lower value of the bootstrap confidence intervals.
}
\item{auc_ci_upper}{Upper value of the bootstrap confidence intervals.
}

If multiple imputation was used the output also contains an \code{mids_analyses} object, containing one \code{adjusted_rmst} object for each imputed dataset.

}
\references{
Sarah C. Conner, Lisa M. Sullivan, Emelia J. Benjamin, Michael P. LaValley, Sandro Galea, and Ludovic Trinquart (2019). "Adjusted Restricted Mean Survival Times in Observational Studies". In: Statistics in Medicine 38, pp. 3832-3860

Patrick Royston and Mahesh K. B. Parmar (2013). "Restricted Mean Survival Time: An Alternative to the Hazard Ratio for the Design and Analysis of Randomized Trials with a Time-To-Event Outcome". In: BMC Medical Research Methodology 13.152
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
# including standard errors and confidence intervals
adjrmst <- adjusted_rmst(adjsurv, from=0, to=1, use_boot=TRUE,
                         conf_level=0.95)
print(adjrmst)
}
}
