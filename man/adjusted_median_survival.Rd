\name{adjusted_median_survival}
\alias{adjusted_median_survival}

\title{
Calculate Confounder-Adjusted Median Survival Times
}
\description{
This function can be utilized to calculate the confounder-adjusted median survival time, given previously estimated adjusted survival curves.
}
\usage{
adjusted_median_survival(adjsurv, use_boot=F, conf_level=0.95, verbose=T)
}

\arguments{
  \item{adjsurv}{
An \code{adjustedsurv} object created using the \code{adjustedsurv} function.
  }
  \item{use_boot}{
Whether bootstrap estimates should be used to calculate the confidence intervals of the estimates. Can only be used if \code{bootstrap=TRUE} was used in the \code{\link{adjustedsurv}} call.
  }
  \item{conf_level}{
A number specifying the confidence level of the bootstrap confidence intervals.
  }
  \item{verbose}{
Whether to print the output or not.
  }
}
\details{

The median survival time is simply the time when half the patients are expected to be alive. The chance of surviving beyond that time is 50 percent. It can be read directly from the respective survival curve by drawing a straight line from the 0.5 point on the Y-Axis and reading the Y-Axis value where this line intersects with the survival curve. If the survival probability never drops below 50 percent, the median survival time cannot be calculated.

This function calculates this quantity automatically. Confidence intervals can be estimated using bootstrapping.

}
\value{

Silently returns the \code{data.frame} object which can be seen when \code{verbose=TRUE}. The column \code{median_surv} includes the confounder-adjusted median survival times, the columns \code{ci_lower} and \code{ci_upper} include the lower and upper end of the bootstrap confidence intervals, if estimated. The column \code{n} shows the number of times the median survival time could be calculated. This might be smaller than \code{n_boot} in the original \code{adjustedsurv} call. In this case the estimates should be treated with caution.

}

\author{
Robin Denz
}

\seealso{
\code{\link{adjustedsurv}}
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

# calculate adjusted median survival times
adjusted_median_survival(adjsurv, use_boot=FALSE)

# calculate adjusted median survival times with confidence intervals
adjusted_median_survival(adjsurv, use_boot=TRUE, conf_level=0.95)

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
