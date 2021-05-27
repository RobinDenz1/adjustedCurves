\name{cif_direct_pseudo}
\alias{cif_direct_pseudo}

\title{
Direct Adjusted CIFs using Pseudo-Values
}
\description{

This page explains the details of estimating direct adjusted cumulative incidence functions using pseudo-values in a competing risks setting. The function is called internally when using \code{method="direct_pseudo"} in the \code{\link[=adjustedcif]{adjustedcif}} function.

It uses essentially the same methodology as the regular direct adjustment method, with the difference that the outcome mechanism is modeled using pseudo-values.

}
\usage{
cif_direct_pseudo(data, variable, ev_time, event, cause, conf_int,
                  conf_level=0.95, times, outcome_vars, model_type="lm",
                  type_time="factor", spline_df=10)
}

\arguments{
  \item{data}{
A \code{data.frame} object containing the needed time-to-event data in standard format.
}
  \item{variable}{
A character string specifying the variable by which the CIFs should be grouped. Must be a valid columname of \code{data}.
}
  \item{ev_time}{
A character string specifying the variable indicating the time-to-event or time-to-censoring. Must be a valid columname of \code{data}.
}
  \item{event}{
A character string specifying the numeric event indicator. The censoring indicator should be coded as 0 and all other events of interest as 1, 2, etc. Must be a valid columname of \code{data}.
  }
  \item{cause}{
The cause of interest for which the cumulative incidence functions should be calculated. Should be a number that appears in the \code{event} column of \code{data}.
  }
  \item{conf_int}{
A logical variable, indicating whether the asymptotic variances and confidence intervals of the CIFs should be calculated. Not available for all methods.
  }
  \item{conf_level}{
A number specifying the confidence level of asymptotic confidence intervals.
  }
  \item{times}{
A numeric vector of time points at which the CIFs should be estimated.
}
  \item{outcome_vars}{
A character vector of column names specifying variables to be used when modeling the outcome mechanism. See details and examples.
}
  \item{model_type}{
The type of modelling used to estimate the survival probabilities. Can be either \code{"lm"} (fitting one linear regression for each time point with the pseudo values as response) or \code{"geese"} (fitting a single generalised estimating equations model using the \code{geese} function).
}
  \item{type_time}{
A character string specifying how the time should be modeled. Only used when \code{model_type="geese"}, ignored otherwise. See details.
}
  \item{spline_df}{
The number of degrees of freedom used for the natural-spline or B-spline function. Only used when \code{model_type="geese"}, also ignored if \code{type_time="factor"}.
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
\item{\strong{Dependencies:} This method relies on the \code{geepack} and \code{\link[prodlim]{prodlim}} packages.}
}

This method works by executing the following steps: (1) First Pseudo-Values for the cause-specific cumulartive incidence function are estimated for each observation in the dataset and some points in time T. Afterwards (2) a new dataset is created in which every individual observation has multiple rows, one for each point in time of interest. (3) This dataset is used to fit a generalised estimating equations (geese) model, using the Pseudo-Values as independent variable. Next (4) multiple copies of the new dataset are created, one for each possible level of the \code{variable} of interest. (5) The \code{variable} is then set to one level for all observations in each dataset. (5) The \code{geese} model is used to predict the CIF at some points in time T for each observation in all dataset copies. (6) Those estimated probabilities are averaged for each dataset at each point in time, resulting in adjusted CIFs for all levels of the group variable at the specified points in time.

It is essentially the same procedure as desribed in \code{\link[=cif_direct]{"direct"}}. The only difference is that instead of relying on a \code{\link[riskRegression]{CSC}} model, this method uses Pseudo-Values and a geese model.

When estimating the \code{geese} model the \code{ev_time} variable is used as a factor by default. This results in one coefficient being estimated for each unique point in time, which can be very slow computationally if there are a lot of unique points in time and/or the dataset has many rows. In these cases it is recommended to use \code{type_time="bs"} or \code{type_time="ns"}, which results in the \code{ev_time} being modeled using B-Splines or Natural Splines. Simulation studies indicate that there is little difference in the estimates when an appropriately large number of \code{spline_df} is used.
}
\value{
Returns a \code{data.frame} object containing the estimated cause-specific adjusted CIFs for some points in time for each level of \code{variable}.
}
\references{
Per K. Andersen, Elisavet Syriopoulou, and Erik T. Parner. "Causal Inference in Survival Analysis using Pseudo-Observations". In: Statistics in Medicine 36 (2017), pp. 2669-2681.

Per K. Andersen and Maja Pohar Perme. "Pseudo-Observations in Survival Analysis". In: Statistical Methods in Medical Research 19 (2010), pp. 71-99.

Aris Perperoglou, Willi Sauerbrei, Michal Abrahamowicz, and Matthias Schmid. "A Review of Spline Function Procedures in R". In: BMC Medical Research Methodology 19 (2019), pp. 1-16
}
\author{
Robin Denz
}

\seealso{
\code{\link[geepack]{geese}}, \code{\link[prodlim]{jackknife}}, \code{\link[splines]{ns}}, \code{\link[splines]{bs}}, \code{\link[riskRegression]{CSC}}
}
\examples{
\dontrun{
# simulate some data as example
sim_dat <- sim_confounded_crisk(n=100, max_t=5)
sim_dat$group <- as.factor(sim_dat$group)

# calculate adjusted CIFs, with time as factor
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      cause=1,
                      method="direct_pseudo",
                      outcome_vars=c("x1", "x2", "x3", "x4", "x5", "x6"),
                      type_time="factor")
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
                      spline_df=5)

# plot the curves
plot(adjcif)
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
