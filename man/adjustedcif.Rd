\name{adjustedcif}
\alias{adjustedcif}

\title{
Calculate Cause-Specific Confounder-Adjusted Cumulative Incidence Functions
}
\description{
This is one of two main functions of this R-Package. It allows the user to calculate cause-specific confounder-adjusted cumulative incidence functions using a variety of different methods. It is basically a wrapper for each \code{cif_} function with some additional functionality, such as bootstrapping. Some of these methods require additional packages to be installed and depending on the specified method there might be additional required arguments in the function call. Documentation on that can be found in the respective \code{cif_} function. See details and examples.
}
\usage{
adjustedcif(data, variable, ev_time, event, cause, method, conf_int=F,
            conf_level=0.95, times=NULL, bootstrap=F, n_boot=500,
            n_cores=1, na.action=options("na.action")[[1]], ...)
}

\arguments{
  \item{data}{
A \code{data.frame} object containing the needed time-to-event data in standard format. Can also be a \code{mids} object created with the \code{mice} package. See details for how this works.
  }
  \item{variable}{
A character string specifying the variable by which the cumulative incidence functions should be grouped. Must be a valid columname of \code{data}. The variable specified should ideally be a factor variable, but integers also work for some methods.
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
  \item{method}{
A character string specifying the adjustment method to use. See details.
  }
  \item{conf_int}{
A logical variable, indicating whether the asymptotic variances and confidence intervals of the cumulative incidence should be calculated. Not available for all methods. More information can be found in the documentation of each method. For an alternative way to get confidence intervals, see the \code{bootstrap} argument.
  }
  \item{conf_level}{
A number specifying the confidence level of asymptotic and/or bootstrap confidence intervals.
  }
  \item{times}{
A numeric vector of time points at which the cumulative incidences should be estimated or \code{NULL}. If \code{NULL} the cumulative incidence is estimated at all points in time at which any event occured in the pooled sample.
  }
  \item{bootstrap}{
A logical variable indicating whether bootstrapping should be performed or not. In bootstrapping, a number of simple random samples with replacement of size \code{nrow(data)} is drawn from \code{data}. For each sample the calculations are repeated and used to estimate standard errors and confidence intervals. This can be used to obtain confidence intervals when asymptotic variance calculations are not available.
  }
  \item{n_boot}{
Number of bootstrap replications to perform. Ignored if \code{bootstrap} is \code{FALSE}.
  }
  \item{n_cores}{
The number of cores to use when calculating bootstrap estimates. Ignored if \code{bootstrap=FALSE}. Is set to 1 by default, resulting in single threaded processing. Internally uses the \code{doParallel} package if \code{n_cores > 1}. In that case it also uses the \code{doRNG} package to make the results replicable. See \code{?doRNG} and \code{?doParallel} for more details. Using multiple cores will speed up the calculation considerably in most cases.
  }
  \item{na.action}{
How missing values should be handled. Can be one of: na.fail, na.omit, na.pass or na.exclude. Also accepts strings of the function names. See \code{?na.action} for more details. By default it uses the na.action which is set in the global options by the respective user.
  }
  \item{...}{
Further arguments passed to the respective \code{cif_method}. See details.
  }
}
\details{
The primary purpose of the \code{adjustedcif} function is to provide a convenient way to calculate confounder-adjusted cumulative incidence functions using any of the methods provided in the literature. A \code{\link[=plot.adjustedcif]{plot}} method is provided to graphically display the estimated cumulative incidence functions as well. Currently the following methods can be used:

\itemize{
  \item{"\link[=cif_direct]{direct}": Direct Standardization based on a Cause-Specific Cox-Regression model.}
  \item{"\link[=cif_direct_pseudo]{direct_pseudo}": Direct Standardization based on Pseudo-Values.}
  \item{"\link[=cif_iptw]{iptw}": A weighted Nelson-Aalen estimator.}
  \item{"\link[=cif_iptw_pseudo]{iptw_pseudo}": A weighted estimator based on Pseudo-Values.}
  \item{"\link[=cif_matching]{matching}": Using propensity score matching to estimate the adjusted CIFs.}
  \item{"\link[=cif_aiptw]{aiptw}": An Augmented Inverse Probability of Treatment Weighting estimator.}
  \item{"\link[=cif_aiptw_pseudo]{aiptw_pseudo}": An Augmented Inverse Probability of Treatment Weighting estimator using Pseudo-Values.}
  \item{"\link[=cif_tmle]{tmle}": Targeted Maximum Likelihood Estimation for CIFs.}
  \item{"\link[=cif_aalen_johansen]{aalen_johansen}": A simple stratified Nelson-Aalen estimator without any form of adjustment.}
}

A short description of each method is contained in the documentation of the respective \code{cif_} function. For more detailed descriptions the cited literature in that same documentation can be used. The documentation for \code{method="direct"} for example can be accessed using \code{?cif_direct}.

For most methods approximations for the asymptotic variance of point estimates of the CIF have been proposed in the literature. Where available, those can be calculated and added to the output object using \code{conf_int=TRUE}. It is however recommended to use bootstrapping to estimate the variance instead, which can be done by setting \code{bootstrap=TRUE}. The \code{n_boot} argument is set to 500 by default. This number was choosen because it worked well in simulations but it does not guarantee convergence in practice. Users are recommended to inspect the bootstrapped estimates and adjust the number of replications accordingly. To allow faster bootstrapping the user can choose to run the function on multiple CPU cores in parallel using the \code{n_cores} argument.

This function also allows multiple imputation via the \code{mice} package. Instead of supplying a single data.frame, the user should create a \code{mids} object using the \code{mice} function and directly pass this to the \code{data} argument. When methods are used which rely on previously estimated treatment or outcome models such as \code{direct} or \code{aiptw}, the user is required to supply a \code{mira} object instead of a single model. In other words: the models have to be fit on every imputed dataset before supplying them to this function. See \code{?mice} and the associated documentation for more information on how to use multiple imputation. When using \code{bootstrap=TRUE} and multiple imputation, the bootstrapping is performed on every imputed dataset separatly. Cumulative Incidences are simply averaged across the imputed datasets according to Rubins Rule. Confidence intervals are calculated by first averaging the standard errors over all imputed datasets and afterwards using this pooled value to obtain a new confidence interval with the normal approximation.

}
\value{
Returns an \code{adjustedcif} object containing the following objects:

\item{adjcif}{A \code{data.frame} of estimated cumulative incidences for \code{cause} at some points in time for each level of \code{variable}. Depending on the arguments used also includes standard errors and confidence intervals at points in time.
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
\item{boot_data_same_t}{The value of the adjusted CIFs at the same points in time as the original CIFs for each bootstrap sample.
}
\item{boot_adjcif}{The mean CIFs of all bootstrap samples and corresponding standard deviations and percentile confidence intervals.
}

When multiple imputation was used, the function additionally contains a \code{mids_analyses} object, containing the \code{adjustedcif} objects for each imputed dataset.

Some method specific objects might also be contained in the output.

}

\references{
Hier meine eigenen Paper
}
\author{
The function itself was written by Robin Denz, but some \code{cif_method} functions include wrappers for functions written by other people. More information can be found in the respective \code{cif_method} documentation.
}

\seealso{
\code{\link{plot.adjustedcif}}, \code{\link{cif_aiptw}}, \code{\link{cif_aiptw_pseudo}}, \code{\link{cif_direct}}, \code{\link{cif_direct_pseudo}},
\code{\link{cif_iptw}}, \code{\link{cif_iptw_pseudo}}, \code{\link{cif_matching}}, \code{\link{cif_tmle}}, \code{\link{cif_aalen_johansen}}
}
\examples{
\dontrun{
library(survival)
library(riskRegression)

# simulate some example data
sim_dat <- sim_confounded_crisk(n=500, max_t=5)
sim_dat$group <- as.factor(sim_dat$group)

# treatment assignment model
glm_mod <- glm(group ~ x2 + x3 + x5 + x6, data=sim_dat, family="binomial")

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
#       greater than 100. This is only shown here for convenience.
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      cause=1,
                      method="direct",
                      outcome_model=cox_mod,
                      conf_int=TRUE,
                      bootstrap=TRUE,
                      n_boot=100)

# using bootstrapping with direct adjustment, run in parallel
# all but one availble processor cores
library(foreach)
library(parallel)
library(doRNG)

adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      cause=1,
                      method="direct",
                      outcome_model=cox_mod,
                      conf_int=TRUE,
                      bootstrap=TRUE,
                      n_boot=100,
                      n_cores=parallel::detectCores()-1)

# More specific examples can be found in the documentation of each
# respective cif_method. See ?cif_ + "method" for more information.
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
