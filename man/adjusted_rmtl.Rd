\name{adjusted_rmtl}
\alias{adjusted_rmtl}

\title{
Calculate Confounder-Adjusted Restricted Mean Time Lost
}
\description{
This function can be utilized to calculate the confounder-adjusted restricted mean time lost, given previously estimated adjusted survival curves / CIFs.
}
\usage{
adjusted_rmtl(adj, to, from=0, conf_int=FALSE,
              conf_level=0.95, interpolation="steps",
              subdivisions=1000)
}

\arguments{
  \item{adj}{
An \code{adjustedsurv} object created using the \code{adjustedsurv} function or a \code{adjustedcif} object created using the \code{adjustedcif} function.
  }
  \item{from}{
A number specifying the left side of the time interval of interest. See details. Usually this should be kept at 0 (default) unless there are good reasons to change it.
  }
  \item{to}{
A number specifying the right side of the time interval of interest. See details.
  }
  \item{conf_int}{
Whether bootstrap estimates should be used to calculate the standard errors and confidence intervals of the RMST estimates. Can only be used if \code{bootstrap=TRUE} was used in the \code{\link{adjustedsurv}} or \code{\link{adjustedcif}} call.
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

The adjusted restricted mean time lost is very similar to the adjusted restricted mean survival time. The only difference is that instead of using the integral of the survival curve, the integral of the cumulative incidence function is taken. This is mostly usefull when there are competing risks and the survival function is not identifiable. It can however also be calculated using single-event survival data. When supplying \code{adjustedsurv} objects the survival curves are simply transformed to CIFs before the integrals are calculated.

This quantity can be interpreted as the mean time it takes to succumb to the event of interest. More information on the method itself can be found in the references. More information on the implementation (concerning missing data, confidence interval calculation, etc.) can be found in the documentation of the \code{\link{adjusted_rmst}} function.

}
\value{
Returns a \code{data.frame} containing the columns \code{group} (groups in \code{variable}) and \code{rmtl} (the estimated restricted mean time lost).

If \code{conf_int=TRUE} was used it additionally contains the columns \code{se} (the standard error of the restricted mean time lost), \code{ci_lower} (lower limit of the confidence interval), \code{ci_upper} (upper limit of the confidence interval) and \code{n_boot} (the actual number of bootstrap estimates used).
}
\references{
Sarah C. Conner and Ludovic Trunquart (2021). "Estimation and Modeling of the Restricted Mean Time Lost in the Presence of Competing Risks". In: Statistics in Medicine
}
\author{
Robin Denz
}

\seealso{
\code{\link{adjustedcif}}, \code{\link{adjustedsurv}}
}
\examples{
\dontrun{
library(survival)

###### when using single-event survival data

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

# calculate adjusted restricted mean time lost estimates from 0 to 1,
# including standard errors and confidence intervals
adjrmst <- adjusted_rmst(adjsurv, from=0, to=1, conf_int=TRUE,
                         conf_level=0.95)

###### when using data with competing-risks

library(riskRegression)
library(prodlim)

# simulate some data as example
sim_dat <- sim_confounded_crisk(n=500)
sim_dat$group <- as.factor(sim_dat$group)

# estimate a cause-specific cox-regression model for the outcome
csc_mod <- CSC(Hist(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6 + group,
               data=sim_dat)

# calculate confounder-adjusted cause-specific CIFs for cause = 1
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="direct",
                      outcome_model=csc_mod,
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=500,
                      cause=1)

# calculate adjusted restricted mean time lost estimates from 0 to 1
# including standard errors and confidence intervals
adjrmtl <- adjusted_rmtl(adjcif, from=0, to=1, conf_int=TRUE)

}
}
