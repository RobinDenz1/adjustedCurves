\name{cif_direct_pseudo}
\alias{cif_direct_pseudo}

\title{
Direct Adjusted CIFs using Pseudo-Values
}
\description{

This page explains the details of estimating direct adjusted cumulative incidence functions using pseudo-values in a competing risks setting (\code{method="direct_pseudo"} in the \code{\link[=adjustedcif]{adjustedcif}} function). All regular arguments of the \code{adjustedcif} function can be used. Additionally, the \code{outcome_vars} argument has to be specified in the \code{adjustedcif} call. Further arguments specific to this method are listed below.

}

\arguments{
  \item{outcome_vars}{
[\strong{required}] A character vector of column names specifying variables to be used when modeling the outcome mechanism. See details and examples.
}
  \item{type_time}{
A character string specifying how the time should be modeled. Possible values are \code{"factor"} (modeling each point in time as a separate variable, the default), \code{"bs"} (modeling time using B-Splines) or \code{"ns"} (modeling time using natural splines).
}
  \item{spline_df}{
The number of degrees of freedom used for the natural-spline or B-spline function. Defaults to 5. Ignored if \code{type_time="factor"}.
}
}
\details{

\itemize{
\item{\strong{Type of Adjustment:} Requires a character vector of variable names used to model the outcome mechanism (internally uses \code{\link[geepack]{geese}}).}
\item{\strong{Doubly-Robust:} Estimates are not Doubly-Robust.}
\item{\strong{Categorical groups:} Any number of levels in \code{variable} are allowed. Must be a factor variable.}
\item{\strong{Approximate Variance:} Calculations to approximate the variance and confidence intervals are not available. Bootstrapping can still be used to estimate the confidence intervals (see \code{?adjustedcif}).}
\item{\strong{Allowed Time Values:} Allows both continuous and integer time.}
\item{\strong{Bounded Estimates:} Estimates are not guaranteed to be bounded in the 0 to 1 probability range.}
\item{\strong{Monotone Function:} Estimates are not guaranteed to be monotone.}
\item{\strong{Dependencies:} This method relies on the \pkg{geepack} and \pkg{prodlim} packages.}
}

This method works by executing the following steps: (1) First Pseudo-Values for the cause-specific cumulative incidence function are estimated for each observation in the dataset and some points in time T. Afterwards (2) a new dataset is created in which every individual observation has multiple rows, one for each point in time of interest. (3) This dataset is used to fit a generalized estimating equations (geese) model, using the Pseudo-Values as independent variable. Next (4) multiple copies of the new dataset are created, one for each possible level of the \code{variable} of interest. (5) The \code{variable} is then set to one level for all observations in each dataset. (5) The \code{\link[geepack]{geese}} model is used to predict the CIF at some points in time T for each observation in all dataset copies. (6) Those estimated probabilities are averaged for each dataset at each point in time, resulting in adjusted CIFs for all levels of the group variable at the specified points in time.

It is essentially the same procedure as described in \code{\link[=cif_direct]{"direct"}}. The only difference is that instead of relying on a \code{\link[riskRegression]{CSC}} model, this method uses Pseudo-Values and a geese model.

When estimating the \code{geese} model the \code{ev_time} variable is used as a factor by default. This results in one coefficient being estimated for each unique point in time, which can be very slow computationally if there are a lot of unique points in time and/or the dataset has many rows. In these cases it is recommended to use \code{type_time="bs"} or \code{type_time="ns"}, which results in the \code{ev_time} being modeled using B-Splines or Natural Splines. Simulation studies indicate that there is little difference in the estimates when an appropriately large number of \code{spline_df} is used.
}
\value{
Adds the following additional objects to the output of the \code{adjustedcif} function:

\itemize{
\item \code{pseudo_values}: The matrix of estimated pseudo-values.
\item \code{geese_model}: The geese model used to make the predictions.
}
}
\references{
Per Kragh Andersen, Elisavet Syriopoulou, and Erik T. Parner (2017). "Causal Inference in Survival Analysis using Pseudo-Observations". In: Statistics in Medicine 36, pp. 2669-2681

Per Kragh Andersen and Maja Pohar Perme (2010). "Pseudo-Observations in Survival Analysis". In: Statistical Methods in Medical Research 19, pp. 71-99

Aris Perperoglou, Willi Sauerbrei, Michal Abrahamowicz, and Matthias Schmid (2019). "A Review of Spline Function Procedures in R". in: BMC Medical Research Methodology 19.46, pp. 1-16
}
\author{
Robin Denz
}

\seealso{
\code{\link[geepack]{geese}}, \code{\link[prodlim]{jackknife}}, \code{\link[splines]{ns}}, \code{\link[splines]{bs}}, \code{\link[riskRegression]{CSC}}
}
\examples{
library(adjustedCurves)

set.seed(42)

# simulate some data as example
sim_dat <- sim_confounded_crisk(n=30, max_t=1.3)
sim_dat$group <- as.factor(sim_dat$group)

# calculate adjusted CIFs, with time as factor
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      cause=1,
                      method="direct_pseudo",
                      outcome_vars=c("x1", "x2", "x3", "x4", "x5", "x6"),
                      type_time="factor",
                      force_bounds=TRUE,
                      iso_reg=TRUE)
plot(adjcif)

# with time modelled as B-Spline using 5 degrees of freedom
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      cause=1,
                      method="direct_pseudo",
                      outcome_vars=c("x1", "x2", "x3", "x4", "x5", "x6"),
                      type_time="bs",
                      spline_df=5,
                      force_bounds=TRUE,
                      iso_reg=TRUE)

# plot the curves
plot(adjcif)
}
