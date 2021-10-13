\name{surv_direct_pseudo}
\alias{surv_direct_pseudo}

\title{
Direct Adjusted Survival Curves using Pseudo-Values
}
\description{

This page explains the details of estimating direct adjusted survival curves using pseudo-values for single event time-to-event data. The function is called internally when using \code{method="direct_pseudo"} in the \code{\link[=adjustedsurv]{adjustedsurv}} function. It is highly recommended to use \code{\link[=adjustedsurv]{adjustedsurv}} instead of calling this function directly.

}
\usage{
surv_direct_pseudo(data, variable, ev_time, event, conf_int=FALSE,
                   conf_level=0.95, times, outcome_vars,
                   type_time="factor", spline_df=5,
                   censoring_vars=NULL, ipcw_method="binder")
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
Currently not used because asymptotic variance calculations for this estimator have not been proposed yet. Bootstrapping can still be used to estimate the confidence intervals (see \code{?adjustedsurv}).
  }
  \item{conf_level}{
A number specifying the confidence level of asymptotic confidence intervals.
  }
  \item{times}{
A numeric vector of time points at which the survival probability should be estimated.
}
  \item{outcome_vars}{
A character vector of column names specifying variables to be used when modeling the outcome mechanism. See details and examples.
}
  \item{type_time}{
A character string specifying how the time should be modeled.
}
  \item{spline_df}{
The number of degrees of freedom used for the natural-spline or B-spline function. Ignored if \code{type_time="factor"}.
}
  \item{censoring_vars}{
An optional character vector specifying variables in \code{data}. Those are used in the calculation of inverse probability of censoring weighted pseudo observations. See \code{?pseudo_aareg} for more information. Set to \code{NULL} (default) to use standard pseudo-values without corrections for dependent censoring instead.
}
  \item{ipcw_method}{
The specific method used in the calculation of inverse probability of censoring weighted pseudo observations. Can be either \code{"binder"} (default) or \code{"hajek"}. See \code{?pseudo_aareg} for more information. Ignored if \code{censoring_vars=NULL}.
}

}
\details{

\itemize{
\item{\strong{Type of Adjustment:} Requires a character vector of variable names used to model the outcome mechanism (internally uses \code{\link[geepack]{geese}}).}
\item{\strong{Doubly-Robust:} Estimates are not Doubly-Robust.}
\item{\strong{Categorical groups:} Any number of levels in \code{variable} are allowed. Must be a factor variable.}
\item{\strong{Approximate Variance:} Calculations to approximate the variance and confidence intervals are not available. Users need to use bootstrapping if those are desired.}
\item{\strong{Allowed Time Values:} Allows both continuous and integer time.}
\item{\strong{Bounded Estimates:} Estimates are not guaranteed to be bounded in the 0 to 1 probability range.}
\item{\strong{Monotone Function:} Estimates are not guaranteed to be monotone.}
\item{\strong{Dependencies:} This method relies on the \pkg{geepack} and \pkg{prodlim} packages. Additionally requires the \pkg{eventglm} package if \code{censoring_vars} is specified.}
}

This method works by executing the following steps: (1) First Pseudo-Values for the survival probabilities are estimated for each observation in the dataset and some points in time T. Afterwards (2) a new dataset is created in which every individual observation has multiple rows, one for each point in time of interest. (3) This dataset is used to fit a generalised estimating equations (geese) model, using the Pseudo-Values as independent variable. Next (4) multiple copies of the new dataset are created, one for each possible level of the \code{variable} of interest. (5) The \code{variable} is then set to one level for all observations in each dataset. (5) The \code{geese} model is used to predict the survival probabilities at some points in time T for each observation in all dataset copies. (6) Those estimated probabilities are averaged for each dataset at each point in time, resulting in adjusted survival probabilities for all levels of the group variable at the specified points in time.

It is essentially the same procedure as described in \code{\link[=surv_direct]{"direct"}}. The only difference is that instead of relying on a \code{coxph} model, this method uses Pseudo-Values and a geese model. This can be useful if the data does not conform to some assumptions needed to use the Cox-Regression model (for example the proportional hazards assumption).

When estimating the \code{geese} model the \code{ev_time} variable is used as a factor by default. This results in one coefficient being estimated for each unique point in time, which can be very slow computationally if there are a lot of unique points in time and/or the dataset has many rows. In these cases it is recommended to use \code{type_time="bs"} or \code{type_time="ns"}, which results in the \code{ev_time} being modeled using B-Splines or Natural Splines. Simulation studies indicate that there is little difference in the estimates when an appropriately large number of \code{spline_df} is used.

Additionally, covariate-dependent censoring can be corrected for using inverse probability of censoring weighted pseudo-values (Binder et al. 2014) instead of regular pseudo-values (specified using the \code{censoring_vars} and \code{ipcw_method} arguments).
}
\value{
Returns a \code{list} object containing a \code{data.frame} with the estimated adjusted survival probabilities for some points in time for each level of \code{variable}, the estimated pseudo-values and the geese model.
}
\references{
Per K. Andersen, Elisavet Syriopoulou, and Erik T. Parner. "Causal Inference in Survival Analysis using Pseudo-Observations". In: Statistics in Medicine 36 (2017), pp. 2669-2681.

Per K. Andersen and Maja Pohar Perme. "Pseudo-Observations in Survival Analysis". In: Statistical Methods in Medical Research 19 (2010), pp. 71-99.

Aris Perperoglou, Willi Sauerbrei, Michal Abrahamowicz, and Matthias Schmid. "A Review of Spline Function Procedures in R". In: BMC Medical Research Methodology 19 (2019), pp. 1-16

Nadine Binder, Thomas A. Gerds, and Per Kragh Andersen. "Pseudo-Observations for Competing Risks with Covariate Dependent Censoring". In: Lifetime Data Analysis 20 (2014), pp. 303–315.
}
\author{
Robin Denz
}

\seealso{
\code{\link[geepack]{geese}}, \code{\link[prodlim]{jackknife}}, \code{\link[splines]{ns}}, \code{\link[splines]{bs}}
}
\examples{
\dontrun{
# simulate some data as example
sim_dat <- sim_confounded_surv(n=500, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# calculate adjusted survival curves, with time as factor
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="direct_pseudo",
                        outcome_vars=c("x1", "x2", "x3", "x4", "x5", "x6"),
                        type_time="factor")

# with time modelled as B-Spline using 5 degrees of freedom
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="direct_pseudo",
                        outcome_vars=c("x1", "x2", "x3", "x4", "x5", "x6"),
                        type_time="bs",
                        spline_df=5)

# plot the curves
plot(adjsurv)
}
}
