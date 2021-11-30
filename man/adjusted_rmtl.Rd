\name{adjusted_rmtl}
\alias{adjusted_rmtl}

\title{
Calculate Confounder-Adjusted Restricted Mean Time Lost
}
\description{
This function can be utilized to calculate the confounder-adjusted restricted mean time lost, given previously estimated adjusted survival curves / CIFs.
}
\usage{
adjusted_rmtl(adj, to, from=0, use_boot=FALSE,
              conf_level=0.95)
}

\arguments{
  \item{adj}{
An \code{adjustedsurv} object created using the \code{adjustedsurv} function or a \code{adjustedcif} object created using the \code{adjustedcif} function.
  }
  \item{from}{
A number specifying the left side of the time interval of interest. See details.
  }
  \item{to}{
A number specifying the right side of the time interval of interest. See details.
  }
  \item{use_boot}{
Whether bootstrap estimates should be used to calculate the standard deviation and confidence intervals of the RMST estimates. Can only be used if \code{bootstrap=TRUE} was used in the \code{\link{adjustedsurv}} or \code{\link{adjustedcif}} call.
  }
  \item{conf_level}{
A number specifying the confidence level of the bootstrap confidence intervals.
  }
}
\details{

The adjusted restricted mean time lost is very similar to the adjusted restricted mean survival time. The only difference is that instead of using the integral of the survival curve, the integral of the cumulative incidence function is taken. This is mostly usefull when there are competing risks and the survival function is not identifiable. It can however also be calculated using single-event survival data. When supplying \code{adjustedsurv} objects the survival curves are simply transformed to CIFs before the integrals are calculated.

This quantity can be interpreted as the mean time it takes to succumb to the event of interest. More information on the method itself can be found in the references. More information on the implementation (concerning multiple imputation etc.) can be found in the documentation of the \code{\link{adjusted_rmst}} function.

}
\value{
Returns a \code{adjusted_rmtl} object containing the following object:

\item{auc}{A named vector containing the adjusted restricted mean time lost estimates.
}
\item{from}{\code{from} argument in function call.
}
\item{to}{\code{to} argument in function call.
}

If \code{use_boot} is set to \code{TRUE} and \code{bootstrap} was also set to \code{TRUE} in the original \code{adjustedsurv}/\code{adjustedcif} function call, the following objects are also added to the output:

\item{conf_level}{The confidence level used when calculating the confidence intervals.
}
\item{n_boot}{The numbers of bootstrap iterations actually used to calculate the bootstrap estimates.
}
\item{booted_auc}{A \code{data.frame} containing the adjusted restricted mean time lost estimates for each bootstrap sample.
}
\item{auc_se}{A named vector containing the bootstrapped standard errors of the restricted mean time lost.
}
\item{auc_ci_lower}{Lower value of the bootstrap confidence intervals.
}
\item{auc_ci_upper}{Upper value of the bootstrap confidence intervals.
}

If multiple imputation was used the output also contains an \code{mids_analyses} object, containing one \code{adjusted_rmtl} object for each imputed dataset.

}
\references{
Sarah C. Conner and Ludovic Trunquart (2021). "Estimation and Modeling of the Restricted Mean Time Lost in the Presence of Competing Risks". In: Statistics in Medicine
}
\author{
Robin Denz
}

\seealso{
\code{\link{print.adjusted_rmtl}}, \code{\link{adjustedcif}}, \code{\link{adjustedsurv}}
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
adjrmst <- adjusted_rmst(adjsurv, from=0, to=1, use_boot=FALSE)
print(adjrmst)

# calculate adjusted restricted mean time lost estimates from 0 to 1,
# including standard errors and confidence intervals
adjrmst <- adjusted_rmst(adjsurv, from=0, to=1, use_boot=TRUE,
                         conf_level=0.95)
print(adjrmst)

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
adjrmtl <- adjusted_rmtl(adjcif, from=0, to=1, use_boot=TRUE)
print(adjrmtl)

}
}
