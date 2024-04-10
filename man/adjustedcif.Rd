\name{adjustedcif}
\alias{adjustedcif}

\title{
Estimate Cause-Specific Confounder-Adjusted Cumulative Incidence Functions
}
\description{
This is one of two main functions of this R-Package. It allows the user to estimate cause-specific confounder-adjusted cumulative incidence functions in the presence of competing events using a variety of different methods. Some of these methods require additional packages to be installed and, depending on the specified method, there might be additional required arguments in the function call. More information is available on the documentation page of the respective \code{cif_method}.
}
\usage{
adjustedcif(data, variable, ev_time, event, cause,
            method, conf_int=FALSE, conf_level=0.95,
            times=NULL, bootstrap=FALSE, n_boot=500,
            n_cores=1, na.action=options()$na.action,
            clean_data=TRUE, iso_reg=FALSE,
            force_bounds=FALSE, mi_extrapolation=FALSE,
            ...)
}

\arguments{
  \item{data}{
A \code{data.frame} object containing the needed time-to-event data in standard format. Can also be a \code{mids} object created with the \pkg{mice} package. See details for how this works.
  }
  \item{variable}{
A character string specifying the variable by which the cumulative incidence functions should be grouped. Must be a valid column name of \code{data}. The variable specified should needs to be a factor variable.
  }
  \item{ev_time}{
A character string specifying the variable indicating the time-to-event or time-to-censoring. Must be a valid column name of \code{data}.
  }
  \item{event}{
A character string specifying the numeric event indicator. The censoring indicator should be coded as 0 and all other events of interest as 1, 2, etc. Must be a valid column name of \code{data}.
  }
  \item{cause}{
The cause of interest for which the cumulative incidence functions should be estimated. Should be a number that appears in the \code{event} column of \code{data}.
  }
  \item{method}{
A character string specifying the adjustment method to use. Case sensitive. See details.
  }
  \item{conf_int}{
A logical variable, indicating whether the asymptotic variances and confidence intervals of the cumulative incidence should be estimated. Not available for all methods. More information can be found in the documentation of each method. For an alternative way to get confidence intervals, see the \code{bootstrap} argument.
  }
  \item{conf_level}{
A number specifying the confidence level of asymptotic and/or bootstrap confidence intervals.
  }
  \item{times}{
A numeric vector of time points at which the cumulative incidences should be estimated or \code{NULL}. If \code{NULL} the cumulative incidence is estimated at all points in time at which any event occurred in the pooled sample.
  }
  \item{bootstrap}{
A logical variable indicating whether bootstrapping should be performed or not. In bootstrapping, a number of simple random samples with replacement of size \code{nrow(data)} are drawn from \code{data}. For each sample the calculations are repeated and used to estimate standard errors and confidence intervals. This can be used to obtain confidence intervals when asymptotic variance calculations are not available.
  }
  \item{n_boot}{
Number of bootstrap replications to perform. Ignored if \code{bootstrap} is \code{FALSE}.
  }
  \item{n_cores}{
The number of cores to use when calculating bootstrap estimates. Ignored if \code{bootstrap=FALSE}. Is set to 1 by default, resulting in single threaded processing. Internally uses the \pkg{doParallel} package if \code{n_cores > 1}. In that case it also uses the \pkg{doRNG} package to make the results replicable. See \code{?doRNG} and \code{?doParallel} for more details. Using multiple cores will speed up the calculation considerably in most cases.
  }
  \item{na.action}{
How missing values should be handled. Can be one of: na.fail, na.omit, na.pass or na.exclude. Also accepts strings of the function names. See \code{?na.action} for more details. By default it uses the na.action which is set in the global options by the respective user.
  }
  \item{clean_data}{
If \code{TRUE} all columns which are not needed for the estimation are removed from \code{data} before any further calculations are performed. This ensures that calls to \code{na.omit} (see argument \code{na.action}) do not remove rows which are fully observed in respect to relevant columns due to missing values in irrelevant columns. Set to \code{FALSE} to skip this step. Usually this argument can be ignored. When using non-standard outcome models however it should be set to \code{FALSE}.
  }
  \item{iso_reg}{
Either \code{TRUE} or \code{FALSE} (default), controlling whether isotonic regression is performed on the resulting failure probability estimates. This can be used to ensure that the CIFs are non-decreasing. Since only a few methods may exhibit this problem, this argument is only relevant for some methods (see method specific documentation).
  }
  \item{force_bounds}{
Either \code{TRUE} or \code{FALSE} (default), controlling whether the resulting failure probability estimates should be forced to lie between 0 and 1. If \code{TRUE} and there are values higher than 1, they are simply set to 1. Values lower than 0 are similarly set to 0. Since only a few methods may exhibit this problem, this argument is only relevant for some methods (see method specific documentation).
  }
  \item{mi_extrapolation}{
Whether to allow extrapolation due to imputed survival times or not. This argument is only relevant when using multiply imputed \code{data} with missing covariates in \code{variable}, \code{ev_time} or \code{event}. Depending on the algorithm used to obtain the imputed datasets, it may be possible that one or more imputed datasets contain survival times in a group that are larger than the maximum observed survival time in that group. This may lead to unwanted extrapolation (e.g. the survival curves extending further than they should). By keeping this argument at \code{FALSE}, these times are removed from the output. If set to \code{TRUE}, all available estimates will be used.
  }
  \item{...}{
Further arguments passed to the respective \code{cif_method}. For example when using \code{method="direct"} all further arguments are passed to the \code{cif_direct} function. See details.
  }
}
\details{
The primary purpose of the \code{adjustedcif} function is to provide a convenient way to estimate confounder-adjusted cumulative incidence functions using any of the methods provided in the literature. A \code{\link[=plot.adjustedcif]{plot}} method is provided to graphically display the estimated cumulative incidence functions as well. Currently the following methods can be used:

\itemize{
  \item{"\link[=cif_direct]{direct}": Direct Standardization based on a model describing the outcome mechanism (\code{\link[riskRegression]{CSC}}, \code{\link[riskRegression]{FGR}}, ..).}
  \item{"\link[=cif_direct_pseudo]{direct_pseudo}": Direct Standardization based on Pseudo-Values.}
  \item{"\link[=cif_iptw]{iptw}": A weighted Aalen-Johansen estimator.}
  \item{"\link[=cif_iptw_pseudo]{iptw_pseudo}": A weighted estimator based on Pseudo-Values.}
  \item{"\link[=cif_matching]{matching}": Using Propensity Score Matching to estimate the adjusted CIFs.}
  \item{"\link[=cif_aiptw]{aiptw}": An Augmented Inverse Probability of Treatment Weighting estimator.}
  \item{"\link[=cif_aiptw_pseudo]{aiptw_pseudo}": An Augmented Inverse Probability of Treatment Weighting estimator using Pseudo-Values.}
  \item{"\link[=cif_aalen_johansen]{aalen_johansen}": A simple stratified Aalen-Johansen estimator without any form of adjustment.}
}

A short description of each method is contained in the documentation of the respective \code{cif_method} function. A concise overview of the supported functionality of each method can be found in the associated vignette (\code{vignette(topic="method_overview", package="adjustedCurves")}). For more detailed descriptions the cited literature in the respective documentation pages can be used. The documentation for \code{method="direct"} for example can be accessed using \code{?cif_direct}.

\strong{\emph{Required & Optional Arguments}}

Every method requires the specification of the \code{data}, \code{variable}, \code{ev_time}, \code{event}, \code{cause} and \code{method} arguments. All other arguments mentioned on this page are optional and work for all methods. Depending on the method used, other arguments are required as well. Those can be found on the top of the help page of the respective method. The help pages also list additional optional arguments.

\strong{\emph{Confidence Intervals}}

For most methods approximations for the asymptotic variance of point estimates of the CIF have been proposed in the literature. Where available, those can be estimated and added to the output object using \code{conf_int=TRUE}. It is however recommended to use bootstrapping to estimate the variance instead, which can be done by setting \code{bootstrap=TRUE}. The \code{n_boot} argument is set to 500 by default. This number was chosen because it worked well in simulations but it does not guarantee convergence in practice. Users are recommended to inspect the bootstrapped estimates and adjust the number of replications accordingly. To allow faster bootstrapping the user can choose to run the function on multiple CPU cores in parallel using the \code{n_cores} argument.

\strong{\emph{Missing Data}}

There are two ways to deal with missing data using this function. The first is using the \code{na.action} argument. It simply calls the respective \code{na.action} function on the data before doing any further processing. By using \code{na.action="na.omit"} for example, only rows with complete data are kept for the analysis.

Alternatively, this function also supports the use of multiple imputation via the \pkg{mice} package. Instead of supplying a single data.frame, the user should create a \code{mids} object using the \code{mice} function and directly pass this to the \code{data} argument. When methods are used which rely on previously estimated treatment or outcome models such as \code{"direct"} or \code{"aiptw"}, the user is required to supply a \code{mira} object instead of a single model. In other words: the models have to be fit on every imputed dataset before supplying them to this function. See \code{?mice} and the associated documentation for more information on how to use multiple imputation. When using \code{bootstrap=TRUE} and multiple imputation, the bootstrapping is performed on every imputed dataset separately. Cumulative Incidences are simply averaged across the imputed datasets according to Rubins Rule. Confidence intervals are calculated by first averaging the standard errors over all imputed datasets and afterwards using this pooled value to obtain a new confidence interval with the normal approximation.

\strong{\emph{Competing Risks}}

This function is meant to be used for data containing multiple competing risks. If the data does not contain competing-events, it is recommended to use the \code{adjustedsurv} function instead. It does not estimate the CIF directly, but the CIF can be calculated from the survival using CIF = 1 - \eqn{S(t)}. This can be done automatically in the \code{plot.adjustedsurv} function using \code{cif=TRUE}.

\strong{\emph{Graphical Displays}}

A general plot of the estimated adjusted CIFs can be obtained using the associated \code{\link[=plot.adjustedcif]{plot}} method. In addition, a plot of the difference between two estimated adjusted CIFs can be produced using the \code{\link{plot_curve_diff}} function.

}
\value{
Returns an \code{adjustedcif} object containing the following objects:

\item{adj}{A \code{data.frame} of estimated cumulative incidences for \code{cause} at some points in time for each level of \code{variable}. Depending on the arguments used also includes standard errors and confidence intervals.
}
\item{method}{
The method used to adjust the CIFs.
}
\item{categorical}{Whether there are more than 2 groups in \code{variable}.
}
\item{call}{
The original function call.
}

When the argument \code{bootstrap} is set to \code{TRUE} it additionally contains the following objects:

\item{boot_data}{The adjusted CIFs estimated in each bootstrap sample.
}
\item{boot_adj}{The mean CIFs of all bootstrap samples and corresponding standard errors and percentile confidence intervals.
}

When multiple imputation was used, the function additionally contains a \code{mids_analyses} object, containing the \code{adjustedcif} objects for each imputed dataset.

Some method specific objects might also be contained in the output.

}

\references{
Robin Denz, Renate Klaaßen-Mielke, and Nina Timmesfeld (2023). "A Comparison of Different Methods to Adjust Survival Curves for Confounders". In: Statistics in Medicine 42.10, pp. 1461-1479

More relevant literature can be found in the respective \code{cif_method} documentation.
}
\author{
The function itself was written by Robin Denz, but some \code{cif_method} functions include wrappers for functions written by other people. More information can be found in the respective \code{cif_method} documentation.
}

\seealso{
\code{\link{plot.adjustedcif}}, \code{\link{adjusted_rmtl}}, \code{\link{plot_curve_diff}}
}
\examples{
library(adjustedCurves)
library(survival)

set.seed(42)

# simulate some example data
sim_dat <- sim_confounded_crisk(n=50, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# treatment assignment model
glm_mod <- glm(group ~ x2 + x3 + x5 + x6, data=sim_dat, family="binomial")

if (requireNamespace("riskRegression")) {

library(riskRegression)

# outcome model
cox_mod <- CSC(Hist(time, event) ~ x1 + x2 + x4 + x5 + group, data=sim_dat)

# using direct adjustment with asymptotic confidence intervals for cause 1
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      cause=1,
                      method="direct",
                      outcome_model=cox_mod,
                      conf_int=TRUE,
                      bootstrap=FALSE)

# using IPTW with asymptotic confidence intervals for cause 2
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw",
                      cause=2,
                      treatment_model=glm_mod,
                      conf_int=TRUE,
                      bootstrap=FALSE)

# using AIPTW with asymptotic confidence intervals for cause 1
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      cause=1,
                      method="aiptw",
                      outcome_model=cox_mod,
                      treatment_model=glm_mod,
                      conf_int=TRUE,
                      bootstrap=FALSE)

# using direct adjustment at custom points in time
custom_times <- c(0.001, 0.1, 0.2, 0.6, 1.1)
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      cause=1,
                      method="direct",
                      outcome_model=cox_mod,
                      conf_int=TRUE,
                      bootstrap=FALSE,
                      times=custom_times)

# using bootstrapping with direct adjustment
# NOTE: In practice the number of bootstrap replications should be
#       greater than 10. This is only shown here for convenience.
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      cause=1,
                      method="direct",
                      outcome_model=cox_mod,
                      conf_int=TRUE,
                      bootstrap=TRUE,
                      n_boot=10)

}

# not run because those are too slow
\donttest{
# using bootstrapping with direct adjustment, run in parallel
# on two cores
library(foreach)
library(parallel)
library(doRNG)

if (requireNamespace("riskRegression")) {

adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      cause=1,
                      method="direct",
                      outcome_model=cox_mod,
                      conf_int=TRUE,
                      bootstrap=TRUE,
                      n_boot=4,
                      n_cores=2)
}

if (requireNamespace("mice") & requireNamespace("WeightIt")) {

# using multiple imputation
library(mice)
library(WeightIt)

# simulate some data as example
sim_dat <- sim_confounded_crisk(n=100, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# introduce random missingness in x1 as example
# NOTE: This is only done as an example, in reality you would
#       already have missing data, not introduce it yourself.
sim_dat$x1 <- ifelse(runif(n=50) < 0.5, sim_dat$x1, NA)

# perform multiple imputation
mids <- mice::mice(data=sim_dat, method="pmm", m=2, printFlag=FALSE)

# IPTW Pseudo using WeightIt on imputed data, for cause = 1
adj <- adjustedcif(data=mids,
                   variable="group",
                   ev_time="time",
                   event="event",
                   method="iptw_pseudo",
                   cause=1,
                   treatment_model=group ~ x1 + x2 + x5 + x6,
                   weight_method="ps")
plot(adj, force_bounds=TRUE, iso_reg=TRUE)
}

# More specific examples can be found in the documentation of each
# respective cif_method. See ?cif_ + "method" for more information.
}
}
