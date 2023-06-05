\name{adjusted_curve_diff}
\alias{adjusted_curve_diff}

\title{
Calculate the difference between two Confounder-Adjusted Survival Curves or CIFs
}
\description{
Given a previously created \code{adjustedsurv} or \code{adjustedcif} object, calculate the difference between two of the \code{variable} specific curves. Can either calculate the whole difference curve or the difference at specified points in time.
}
\usage{
adjusted_curve_diff(adj, group_1=NULL, group_2=NULL,
                    times=NULL, conf_int=FALSE, conf_level=0.95,
                    use_boot=FALSE, interpolation="steps")
}

\arguments{
  \item{adj}{
An \code{adjustedsurv} object created using the \code{adjustedsurv} function, or a \code{adjustedcif} object created using the \code{adjustedcif} function.
  }
  \item{group_1}{
Optional argument to get a specific difference. This argument takes a single character string specifying one of the levels of the \code{variable} used in the original \code{adjustedsurv} or \code{adjustedcif} function call. This group will be subtracted from. For example if \code{group_1="A"} and \code{group_2="B"} the difference \code{A - B} will be used. If \code{NULL}, the order of the factor levels in the original \code{data} determines the order. If not \code{NULL}, the \code{group_2} argument also needs to be specified.
  }
  \item{group_2}{
Also a single character string specifying one of the levels of \code{variable}. This corresponds to the right side of the difference equation. See argument \code{group_1}.
  }
  \item{times}{
An optional numeric vector of points in time at which the difference should be estimated. If \code{NULL} (default) the differences are estimated for the whole curve.
  }
  \item{conf_int}{
Whether standard errors, confidence intervals and p-values should be calculated. Only possible when either \code{conf_int=TRUE} or \code{bootstap=TRUE} was used in the original function call. P-values are calculated using a one-sample t-test with the null-hypothesis being that the difference is equal to 0.
  }
  \item{conf_level}{
A number specifying the confidence level of the confidence intervals.
  }
  \item{use_boot}{
Whether to use the standard errors estimated using bootstrapping for the confidence interval and p-value calculation. Can only be used if \code{bootstrap=TRUE} was used in the original \code{adjustedsurv} or \code{adjustedcif} function call. Ignored if \code{conf_int=FALSE}.
  }
  \item{interpolation}{
Either \code{"steps"} (default) or \code{"linear"}. This parameter controls how interpolation is performed. If this argument is set to \code{"steps"}, the curves will be treated as step functions. If it is set to \code{"linear"}, the curves wil be treated as if there are straight lines between the point estimates instead. Points that lie between estimated points will be interpolated accordingly. Should usually be kept at \code{"steps"}. See Details.
  }
}
\details{

\strong{\emph{Confidence Intervals & P-Values}}

The standard error of the difference is estimated using the pooled standard error of the two probability estimates, given by:
\deqn{SE_{group_1 - group_2} = \sqrt{SE_{group_1}^2 + SE_{group_2}^2}}

Confidence intervals are then calculated using this pooled standard error and the normal approximation. The P-Values are also obtained using this standard error combined with a two-sided one-sample t-test. The null-hypothesis is that the difference is equal to 0, and the alternative hypothesis is that the difference is not equal to 0. If p-values are calculated for multiple points in time simultaneously, the user should adjust those. See \code{?p.adjust} for more information.

\strong{\emph{Overall Difference Test}}

This function does not perform a test of the overall difference between two functions. To calculate the integral of the difference in a given interval the \code{\link{plot_curve_diff}} function can be used. Additionally, to test whether that integral is equal to zero the \code{\link{adjusted_curve_test}} function can be used.

\strong{\emph{More than Two Groups}}

If more than two groups are present in \code{variable}, all other comparisons except for \code{group_1 - group_2} are ignored. If multiple comparisons are desired, the user needs to call this function multiple times and adjust the \code{group_1} and \code{group_2} arguments accordingly.

\strong{\emph{Graphical Displays}}

There is no directly associated \code{plot} method for this function. However, this function is used internally when calling the \code{\link{plot_curve_diff}} function. In order to get a plot of the difference curve or point estimates, that function can be used.

\strong{\emph{Multiple Imputation}}

This function works exactly the same way for adjusted survival curves or adjusted CIFs estimated using multiple imputation as it does without any missing values. If multiple imputation was used previously, this function simply uses the pooled estimates to calculate the differences.

\strong{\emph{Computational Details}}

When estimating the difference at some point in time at which no direct point estimates are available, this function needs to interpolate the curves. The interpolation method can be controlled using the \code{interpolation} function. In most cases, the estimated curves are step functions and the default (\code{interpolation="steps"}) is therefore appropriate. However, when parametric survival models where used in the estimation process it might be preferable to use linear interpolation instead.

}
\value{
Returns a \code{data.frame} containing the columns \code{time} (the points in time where the difference was estimated) and \code{diff} (the estimated difference).

If \code{conf_int=TRUE} was used in the function call, it additionally contains the columns \code{se} (the estimated standard error of the difference), \code{ci_lower} (lower limit of the confidence interval of the difference), \code{ci_upper} (upper limit of the confidence interval of the difference) and \code{p_value} (the p-value for the test of pointwise difference).
}

\references{
John P. Klein, Brent Logan, Mette Harhoff, and Per Kragh Andersen (2007). "Analyzing Survival Curves at a Fixed Point in Time". In: Statistics in Medicine 26, pp. 4505-4519

Michael Coory, Karen E. Lamb, and Michael Sorich (2014). "Risk-Difference Curves can be used to Communicate Time-Dependent Effects of Adjuvant Therapies for Early Stage Cancer". In: Journal of Clinical Epidemiology 67, pp. 966-972
}

\author{
Robin Denz
}

\seealso{
\code{\link{plot_curve_diff}}, \code{\link{adjustedsurv}}, \code{\link{adjustedcif}}
}
\examples{
library(adjustedCurves)
library(survival)
library(cmprsk)

#### Simple Survival Case with adjusted survival curves ####

# simulate some data as example
set.seed(42)
sim_dat <- sim_confounded_surv(n=30, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# propensity score model
ps_mod <- glm(group ~ x1 + x2 + x4 + x5, data=sim_dat, family="binomial")

# use it to calculate adjusted survival curves with bootstrapping
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="iptw_km",
                        treatment_model=ps_mod,
                        conf_int=TRUE,
                        bootstrap=TRUE,
                        n_boot=10) # n_boot should be much higher in reality

# calculate the whole difference curve
adjdiff <- adjusted_curve_diff(adjsurv)

# only some points in time
adjdiff <- adjusted_curve_diff(adjsurv, times=c(0.2, 0.4))

# with confidence intervals, p-values
adjdiff <- adjusted_curve_diff(adjsurv, times=c(0.2, 0.4), conf_int=TRUE)

# using bootstrapping
adjdiff <- adjusted_curve_diff(adjsurv, times=c(0.2, 0.4), conf_int=TRUE,
                               use_boot=TRUE)

#### Competing Risks Case with adjusted CIFs ####
library(riskRegression)
library(prodlim)

# simulate some data as example
sim_dat <- sim_confounded_crisk(n=41, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# estimate a cause-specific cox-regression for the outcome
csc_mod <- CSC(Hist(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6 + group,
               data=sim_dat)

# use it to calculate adjusted CIFs for cause = 1 with bootstrapping
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="direct",
                      outcome_model=csc_mod,
                      conf_int=TRUE,
                      bootstrap=TRUE,
                      n_boot=10,
                      cause=1)

# calculate the whole difference curve
adjdiff <- adjusted_curve_diff(adjcif)

# with confidence intervals
adjdiff <- adjusted_curve_diff(adjcif, conf_int=TRUE)

# only at specific points in time
adjdiff <- adjusted_curve_diff(adjcif, times=c(0.2, 0.4), conf_int=TRUE)
}
