\name{surv_method_iptw_pseudo}
\alias{surv_method_iptw_pseudo}

\title{
Inverse Probability of Treatment Weighted Survival Estimates using Pseudo-Values
}
\description{
This method works by modelling the treatment assignment mechanism. The individual observations are weighted using inverse probability of treatment weights or other balancing weights (see details) and pseudo-values. If the model is correctly specified, unbiased results will be obtained. The main difference to the \code{surv_method_iptw_cox} and \code{surv_method_iptw_km} functions is the use of pseudo-values. See details for more information.
}
\usage{
surv_method_iptw_pseudo(data, variable, ev_time, event, conf_int,
                        conf_level=0.95, times, treatment_model,
                        weight_method="ps", stabilize=T,
                        se_method="cochrane", censoring_vars=NULL,
                        ipcw_method="binder", ...)
}

\arguments{
  \item{data}{
A \code{data.frame} object containing the needed time-to-event data in standard format.
}
  \item{variable}{
A character string specifying the variable by which the survival curves should be grouped. Must be a valid columname of \code{data}.
}
  \item{ev_time}{
A character string specifying the variable indicating the time-to-event or time-to-censoring. Must be a valid columname of \code{data}.
}
  \item{event}{
A character string specifying the binary event indicator. Must be a valid columname of \code{data}.
}
  \item{conf_int}{
A logical variable, indicating whether the asymptotic variances and confidence intervals of the survival probabilities should be calculated.
  }
  \item{conf_level}{
A number specifying the confidence level of asymptotic and/or bootstrap confidence intervals.
  }
  \item{times}{
A numeric vector of time points at which the survival probability should be estimated or \code{NULL}. If \code{NULL} the survival is estimated at all points in time at which an event occured.
}
  \item{treatment_model}{
Must be either a model object with \code{variable} as response variable, a vector of weights or a formula which can be passed to \code{WeightIt}.
}
  \item{weight_method}{
Method used in \code{WeightIt} function call. Ignored if \code{treatment_model} is not a formula object.
}
  \item{stabilize}{
Whether to stabilize the weights or not. Is set to \code{TRUE} by default. Stabilizing has no effect on the estimates, only on the asymptotic variance calculations and confidence intervals. Without stabilizing, those estimates are biased downward.
  }
  \item{se_method}{
One of \code{"miller"}, \code{"galloway"}, \code{"cochrane"} and \code{"simple"}. Specifies which kind of standard error to calculate. Defaults to \code{"cochrane"}. See details.
}
  \item{censoring_vars}{
An optional character vector specifying variables in \code{data}. Those are used in the calculation of inverse probability of censoring weighted pseudo observations. See \code{?pseudo_aareg} for more information. Set to \code{NULL} (default) to use standard pseudo-values without corrections for dependent censoring instead.
}
  \item{ipcw_method}{
The specific method used in the calculation of inverse probability of censoring weighted pseudo observations. Can be either \code{"binder"} (default) or \code{"hajek"}. See \code{?pseudo_aareg} for more information. Ignored if \code{censoring_vars=NULL}.
}
  \item{...}{
Further arguments passed to \code{WeightIt}.
  }
}
\details{

\itemize{
\item{\strong{Type of Adjustment:} Requires a model describing the treatment assignment mechanism. This must be either a \code{\link{glm}} or \code{\link{multinom}} object.}
\item{\strong{Doubly-Robust:} Estimates are not Doubly-Robust.}
\item{\strong{Categorical groups:} Any number of levels in \code{variable} are allowed. Must be a factor variable.}
\item{\strong{Approximate Variance:} Calculations to approximate the variance and confidence intervals are available.}
\item{\strong{Allowed Time Values:} Allows both continuous and integer time.}
\item{\strong{Bounded Estimates:} Estimates are not guaranteed to be bounded in the 0 to 1 probability range.}
\item{\strong{Monotone Function:} Estimates are not guaranteed to be monotone.}
\item{\strong{Dependencies:} This method relies on the \code{\link[prodlim]{prodlim}} package. The \code{WeightIt} package is also required if \code{treatment_model} is a formula object.}
}

Adjusted survival curves are calculated by first estimating appropriate case-weights for each observation in \code{data}. This can be done using inverse probability of treatment weights using the propensity score (usually estimated using a logistic regression model) or by some other method (see \code{?weightit}). Pseudo-Values are then calculated for every observation in \code{data} at some points in time \eqn{T}. Since Pseudo-Values bypass the problem of censoring, a simple weighted average of the Pseudo-Values can be taken for every \eqn{T}. See Andersen et al. (2017) for more details on this method and Andersen and Perme (2010) for more information on Pseudo-Values in general.

The standard error of this estimator can be approximated by calculation a weighted version of the standard error estimator. Interestingly, no exact method exists in the weighted case. Four approximations are implemented which can be choosen using the \code{se_method} argument. The equations for \code{"miller"}, \code{"galloway"} and \code{"cochrane"} are described and compared in Gatz and Smith (1995). \code{"simple"} is just the standard equation with a weight term added and should only be used with stabilized weightd (\code{stabilize=TRUE}). The resulting estimates are conservative in nature, because they ignore the variation introduced by the pseudo-value calculation (Zeng et al. 2021). It is generally recommended to use bootstrap estimates instead.

}
\value{
Returns a \code{data.frame} object containing the estimated adjusted survival probabilities for some points in time for each level of \code{variable}.
}
\references{
Per K. Andersen, Elisavet Syriopoulou, and Erik T. Parner. "Causal Inference in Survival Analysis using Pseudo-Observations". In: Statistics in Medicine 36 (2017), pp. 2669-2681.

Per K. Andersen and Maja Pohar Perme. "Pseudo-Observations in Survival Analysis". In: Statistical Methods in Medical Research 19 (2010), pp. 71-99.

Donald F. Gatz and Luther Smith. "The Standard Error of a Weighted Mean Concentration - I: Bootstrapping Vs Other Methods". In: Atmospheric Environment 29.11 (1995), pp. 1185-1193.

William G. Cochran. Sampling Techniques. Vol. 3. New York: Wiley, 1977.

J. N. Galloway, G. E. Likens, and M. E. Hawley. "Acid Precipitation: Natural Versus Anthropogenic Components". In: Science 226 (1984), pp. 829-831.

J. M. Miller. A Statistical Evaluation of the U.S. Precipitation Chemistry Network. Precipitation Scavenging (edited by Semonin R. G. and Beadle R. W.) pp. 639-659. Available as CONF 74100 from National Technical Information Service, U.S. Dept. of Commerce, Springfiel, VA. 1977.

ZENG PAPER
}
\author{
Robin Denz
}

\seealso{
\code{\link[WeightIt]{WeightIt}}, \code{\link[prodlim]{prodlim}}
}
\examples{
\dontrun{
# simulate some data as example
sim_dat <- sim_confounded_surv(n=500, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# estimate a treatment assignment model
glm_mod <- glm(group ~ x1 + x3 + x5 + x6, data=sim_dat, family="binomial")

# use it to calculate adjusted survival curves
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="iptw_km",
                        treatment_model=glm_mod)

# Alternatively, use custom weights
# In this example we use weights calculated using the propensity score,
# which is equal to using the glm model directly in the function
ps_score <- glm_mod$fitted.values
weights <- ifelse(sim_dat$group==1, 1/ps_score, 1/(1-ps_score))

adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="iptw_km",
                        treatment_model=weights)

# And a third alternative: use the WeightIt package
# here an example with equal results to the ones above:
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="iptw_km",
                        treatment_model=group ~ x1 + x3 + x5 + x6,
                        weight_method="ps")

# here an example using Optimization-Based Weighting:
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="iptw_km",
                        treatment_model=group ~ x1 + x3 + x5 + x6,
                        weight_method="optweight")
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
