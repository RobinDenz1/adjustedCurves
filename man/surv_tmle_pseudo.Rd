\name{surv_tmle_pseudo}
\alias{surv_tmle_pseudo}

\title{
Adjusted Survival Curves using Targeted Maximum Likelihood Estimation based on Pseudo-Values
}
\description{
Calculate confounder-adjusted survival curves using a standard targeted maximum likelihood estimator for continuous outcomes for right-censored data using pseudo-values.
}
\usage{
surv_tmle_pseudo(data, variable, ev_time, event, conf_int,
                 conf_level=0.95, times, outcome_vars,
                 treatment_vars=NULL, SL.trt=NULL, SL.ftime=NULL,
                 treatment_model=NULL, cv_folds=5, censoring_vars=NULL,
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
A logical variable, indicating whether the asymptotic variances and confidence intervals of the survival probabilities should be calculated. Not available for all methods.
  }
  \item{conf_level}{
A number specifying the confidence level of asymptotic confidence intervals.
  }
  \item{times}{
A numeric vector of time points at which the survival probability should be estimated.
}
  \item{outcome_vars}{
A character vector of column names in \code{data}, specifying the variables which should be used when modeling the outcome mechanism using \code{SuperLearner}. The columns specified by \code{variable} and \code{ev_time} are included by default and therefore do not have to be mentioned in this argument.
}
  \item{treatment_vars}{
An optional chracter vector of column names in \code{data}, specifying the variables which should be used when modeling the treatment assignment using \code{SuperLearner}. Alternatively, a \code{glm} model or a vector of previously estimated propensity scores can be supplied to the \code{treatment_model} argument. If both are defined, \code{treatment_model} is used.
}
  \item{SL.trt}{
A character vector or list specification to be passed to the \code{SL.library} option in the call to \code{SuperLearner} for the outcome regression. See \code{?SuperLearner} for more information on how to specify valid \code{SuperLearner} libraries.
}
  \item{SL.ftime}{
A character vector or list specification to be passed to the \code{SL.library} option in the call to \code{SuperLearner} for the treatment assignment regression. See \code{?SuperLearner} for more information on how to specify valid \code{SuperLearner} libraries.
}
  \item{treatment_model}{
Must be either a \code{glm} model object with \code{variable} as response variable or a numeric vector of propensity scores. Alternatively, the arguments \code{SL.trt} and \code{treatment_vars} can be specified to estimate propensity scores using \code{SuperLearner}. See details and examples.
}
  \item{cv_folds}{
The number of cross-validation folds to be used when estimating the models using \code{SuperLearner}.
}
  \item{censoring_vars}{
An optional character vector specifying variables in \code{data}. Those are used in the calculation of inverse probability of censoring weighted pseudo observations. See \code{?pseudo_aareg} for more information. Set to \code{NULL} (default) to use standard pseudo-values without corrections for dependent censoring instead.
}
  \item{ipcw_method}{
The specific method used in the calculation of inverse probability of censoring weighted pseudo observations. Can be either \code{"binder"} (default) or \code{"hajek"}. See \code{?pseudo_aareg} for more information. Ignored if \code{censoring_vars=NULL}.
}
  \item{...}{
Further arguments passed to \code{tmle}.
  }
}
\details{

\itemize{
\item{\strong{Type of Adjustment:} Adjustments are made based on the treatment assignment mechanism and the outcome mechanism. No models can be supplied for the outcome mechanism. The adjustments are made based on \code{\link[SuperLearner]{SuperLearner}} libraries. For the treatment assignment, either a \code{glm} model or \code{\link[SuperLearner]{SuperLearner}} can be used.}
\item{\strong{Doubly-Robust:} Estimates are Doubly-Robust.}
\item{\strong{Categorical groups:} Currently only two groups in \code{variable} are allowed. Must be an integer variable only containing 0 and 1.}
\item{\strong{Approximate Variance:} Calculations to approximate the variance and confidence intervals are available.}
\item{\strong{Allowed Time Values:} In contrast to the other TMLE based methods, this method allows both integer time and continuous times.}
\item{\strong{Bounded Estimates:} Estimates are not guaranteed to be bounded in the 0 to 1 probability range.}
\item{\strong{Monotone Function:} Estimates are not guaranteed to be monotone.}
\item{\strong{Dependencies:} This method relies on the \code{tmle}, \code{SuperLearner} and \code{prodlim} packages. Additionally requires the \code{eventglm} package if \code{censoring_vars} is specified.}
}

As in \link[=surv_direct_pseudo]{direct_pseudo}, \link[=surv_iptw_pseudo]{iptw_pseudo} and \link[=surv_aiptw_pseudo]{aiptw_pseudo} this method utilizes pseudo-values to enable the use of the methodology originally developed for mean outcomes. First, pseudo values are estimated for each person for a range of points in time. These pseudo values are subsequently used to estimate an outcome model with an ensemble of machine learning and other methods (using the \code{SuperLearner} package). Using an additional estimate of the propensity score, the targeted maximum likelihood estimation methodology for estimating arithmetic means is applied directly to the pseudo-values. This makes it possible to directly use continuous survival times with targeted maximum likelihood estimation. MORE FROM MY PAPER.

Additionally, covariate-dependent censoring can be corrected for using inverse probability of censoring weighted pseudo-values (Binder et al. 2014) instead of regular pseudo-values (specified using the \code{censoring_vars} and \code{ipcw_method} arguments).
}
\value{
Returns a \code{data.frame} object containing the estimated adjusted survival probabilities for some points in time for each level of \code{variable}.
}
\references{
# TODO: References

Nadine Binder, Thomas A. Gerds, and Per Kragh Andersen. "Pseudo-Observations for Competing Risks with Covariate Dependent Censoring". In: Lifetime Data Analysis 20 (2014), pp. 303–315.
}
\author{
Robin Denz
}

\seealso{
\code{\link[tmle]{tmle}}, \code{\link[SuperLearner]{SuperLearner}}
}
\examples{
\dontrun{
library(survival)

# simulate some data as example
sim_dat <- sim_confounded_surv(n=500, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# TODO: Add examples
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
