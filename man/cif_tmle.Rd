\name{cif_tmle}
\alias{cif_tmle}

\title{
Targeted Maximum Likelihood Estimation for Continuous Time Competing Events Data
}
\description{

This page explains the details of estimating causal cause-specific cumulative incidence functions in a competing risks setting with targeted maximum likelihood estimation (\code{method="tmle"} in the \code{\link[=adjustedcif]{adjustedcif}} function). All regular arguments of the \code{adjustedcif} function can be used. Additionally, the \code{outcome_model} argument and the \code{treatment_model} argument have to be specified in the \code{adjustedcif} call. Further arguments specific to this method are listed below.

}

\arguments{
  \item{outcome_model}{
[\strong{required}] Should be a list containing at least one Cox model formula for the desired cause. In those Cox model formula, the time variable is always called \code{time} and the status variable should always be called \code{status}. For example, to use just one Cox model for \code{cause = 1} including all variables of the dataset as independent variables the user can use \code{list(Surv(time, status==1) ~ .)}. If \code{cause = 2} is of interest it should be \code{list(Surv(time, status==2) ~ .)} instead. All Cox models defined in that list are fitted to the data and the ensemble of these models is used to provide the initial predictions for the conditional hazard. See details and the documentation of the \pkg{concrete} package for more information.
}
  \item{treatment_model}{
[\strong{required}] A character vector specifying which \pkg{SuperLearner} libraries should be used to obtain an estimate of the propensity score. For example, \code{c("SL.glm", "SL.glmnet")} could be used. See \code{?SuperLearner} for more details.
}
  \item{censoring_model}{
Either \code{NULL} (default) to make no adjustments for dependent censoring, or a list of Cox models as described in the \code{outcome_model} argument. The only difference between this and the \code{outcome_model} argument is that \code{status==0} should be used in the Cox formulas. See below or \code{?formatArguments} in the \pkg{concrete} package for more details.
  }
  \item{cv_args}{
A list of arguments specifying how exactly cross-validation should be performed. Internally passed to the \code{CVArg} argument of the \code{formatArguments} function in the \pkg{concrete} package.
  }
  \item{max_update_iter}{
A single positive integer specifying the maximum iterations performed to obtain the estimates. Defaults to 500. Internally passed to the \code{MaxUpdateIter} argument of the \code{formatArguments} function in the \pkg{concrete} package.
  }
  \item{one_step_eps}{
A single positive number specifying the step size of the tmle updates. Defaults to 0.1. Internally passed to the \code{OneStepEps} argument of the \code{formatArguments} function in the \pkg{concrete} package.
  }
  \item{min_nuisance}{
A single number between 0 and 1 used for truncating the g-related denominator of the clever covariate. Defaults to \code{5/sqrt(nrow(data))/log(nrow(data))}. Internally passed to the \code{MinNuisance} argument of the \code{formatArguments} function in the \pkg{concrete} package.
  }
  \item{verbose}{
Whether to print estimation information of the \code{doConcrete} function in the \pkg{concrete} package. Defaults to \code{FALSE}.
  }
  \item{return_models}{
Whether to add the estimated models for the outcome, treatment, and censoring mechanism to the output object. Defaults to \code{TRUE}.
  }
}
\details{

\itemize{
\item{\strong{Type of Adjustment:} Requires a model describing the treatment assignment mechanism and the outcome mechanism, also allows a model for the censoring mechanism. See details and the \pkg{concrete} package.}
\item{\strong{Doubly-Robust:} Estimates are Doubly-Robust.}
\item{\strong{Categorical groups:} This function currently only allows two levels in \code{variable}.}
\item{\strong{Approximate Variance:} Calculations to approximate the variance and confidence intervals are available.}
\item{\strong{Allowed Time Values:} Allows both continuous and integer time.}
\item{\strong{Bounded Estimates:} Estimates are guaranteed to be bounded in the 0 to 1 probability range.}
\item{\strong{Monotone Function:} Estimates are guaranteed to be monotone.}
\item{\strong{Dependencies:} This method relies on the \pkg{concrete} package, the \pkg{data.table} package and all of their respective dependencies.}
}

\strong{\emph{What it does:}}

This function implements Targeted Maximum Likelihood Estimation (TMLE) for continuously distributed competing-events data as described in Rytgaard et al. (2023) and Rytgaard and van der Laan (2022). The TMLE method is similar to the AIPTW methods included in this package. It also relies on both an outcome model and a treatment model (with an additional optional censoring model) to obtain the counterfactual failure probability estimates. In contrast to the AIPTW methods, however, the estimator uses an iterative approach to obtain the estimates where each update targets the entire cumulative incidence function. As a consequence, the resulting estimates are guaranteed to lie in the 0/1 probability bounds and are also guaranteed to be non-decreasing over time. Simulation studies and theoretical results indicate a good performance of this method in terms of bias and standard errors. See the cited literature for more detailed and more rigorous explanations of the method.

Instead of relying on a single model to obtain the propensity score or the initial conditional hazards estimates, this estimator relies on the \pkg{SuperLearner} framework in conjunction with cross-validation to do this. How cross-validation should be performed may be controlled with the \code{cv_args} argument. The resulting models can be inspected from the output object if \code{return_models} is set to \code{TRUE}.

\strong{\emph{The Implementation:}}

Internally, this function simply calls multiple functions of the \pkg{concrete} package in correct order with appropriate arguments. This wrapper function is limited in the sense that it does not allow dynamic interventions or time-varying variables, which are supported by the \pkg{concrete} package. It is recommended to use the \pkg{concrete} package directly when the user wants to use these features or other specific settings are required.

\strong{\emph{Speed Considerations:}}

This method is very computationally expensive. For medium to large datasets and when considering many different points in time, it will usually take a very long time to execute. If speed is important, we recommend using other methods. Alternatively, user may adjust the \code{times} arguments to target fewer points in time.

}
\note{
A previous version of this package (<= 0.9.1) included a function with the same name, which was removed in version 0.10.0. The old version implemented a TMLE estimator that was only applicable to discrete-time survival data based on the \pkg{survtmle} package, which was removed from CRAN. The new version implements a different estimator. Code using this method for version <= 0.9.1 does NOT work with versions 0.10.2 or higher.
}
\value{
Adds the following additional objects to the output of the \code{adjustedcif} function:

\itemize{
\item \code{concrete_object}: The object returned by the \code{doConcrete} function.
}
}
\references{
Helene C. W. Rytgaard and Mark J. van der Laan (2023). "Targeted Maximum Likelihood Estimation for Causal Inference in Survival and Competing Risks Analysis". In: Lifetime Data Analysis

Helene C. W. Rytgaard and Mark J. van der Laan (2023). "One-Step Targeted Maximum Likelihood Estimation for Targeting Cause-Specific Absolute Risks and Survival Curves". In: Biometrika

Helene C. W. Rytgaard, Frank Eriksson and Mark J. van der Laan (2023). "Estimation of Time-Specific Intervention Effects on Continuously Distributed Time-To-Event Outcomes by Targeted Maximum Likelihood Estimation". In: Biometrics

David Chen, Helene C. W. Rytgaard and Edwin Fong and Jens M. Tarp and Maya L. Petersen and Mark J. van der Laan and Thomas A. Gerds (2023). "concrete: An R Package for Continuous-Time, Competing Risks Targeted Maximum Likelihood Estimation". Available at <https://github.com/imbroglio-dc/concrete> or on CRAN
}
\author{
The wrapper function was written by Robin Denz, but the real estimation functions are all contained in the \pkg{concrete} package, which was written by David Chen. See \code{?doConcrete} for more information.
}

\seealso{
\code{\link{adjustedcif}}
}
\examples{
library(adjustedCurves)

data <- sim_confounded_crisk(n=100)
data$group <- factor(data$group)

# for a single point in time using only one model for both
# the treatment mechanism and outcome mechanism
out <- adjustedcif(data=data,
                   variable="group",
                   ev_time="time",
                   event="event",
                   cause=1,
                   treatment_model=c("SL.glm"),
                   outcome_model=list(Surv(time, status==1) ~ .),
                   times=c(0.5),
                   conf_int=TRUE,
                   method="tmle")

## using multiple models for both the treatment assignment and
## outcome mechanism
out <- adjustedcif(data=data,
                   variable="group",
                   ev_time="time",
                   event="event",
                   cause=1,
                   treatment_model=c("SL.glm", "SL.mean"),
                   outcome_model=list(Surv(time, status==1) ~ x1 + x3,
                                      Surv(time, status==1) ~ x2 + x4 + x5),
                   times=c(0.5),
                   conf_int=TRUE,
                   method="tmle")

## with corrections for covariate dependent censoring
out <- adjustedcif(data=data,
                   variable="group",
                   ev_time="time",
                   event="event",
                   cause=1,
                   treatment_model=c("SL.glm", "SL.mean"),
                   outcome_model=list(Surv(time, status==1) ~ x1 + x3,
                                      Surv(time, status==1) ~ x2 + x4 + x5),
                   censoring_model=list(Surv(time, status==0) ~ x6 + x1),
                   times=c(0.5),
                   conf_int=TRUE,
                   method="tmle")
}
