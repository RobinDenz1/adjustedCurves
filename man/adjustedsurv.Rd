\name{adjustedsurv}
\alias{adjustedsurv}

\title{
Calculate Confounder-Adjusted Survival Curves
}
\description{
This is one of the two main functions of this R-Package. It allows the user to calculate confounder-adjusted survival curves using a variety of different methods. It is basically a wrapper for each \code{surv_} function with some additional functionality, such as bootstrapping. Some of these methods require additional packages to be installed and depending on the specified method there might be additional required arguments in the function call. Documentation on that can be found in the respective \code{surv_} function. See details and examples.
}
\usage{
adjustedsurv(data, variable, ev_time, event, method, conf_int=F,
             conf_level=0.95, times=NULL, bootstrap=F, n_boot=500,
             n_cores=1, na.action=options("na.action")[[1]], ...)
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
  \item{method}{
A character string specifying the adjustment method to use. See details.
  }
  \item{conf_int}{
A logical variable, indicating whether the asymptotic variances and confidence intervals of the survival probabilities should be calculated. Not available for all methods.
  }
  \item{conf_level}{
A number specifying the confidence level of asymptotic and/or bootstrap confidence intervals.
  }
  \item{times}{
A numeric vector of time points at which the survival probability should be estimated or \code{NULL} (default). If \code{NULL} the survival is estimated at all points in time at which an event occured in the pooled sample.
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
How missing values should be handled. Can be one of: na.fail, na.omit, na.pass or na.exclude. Also accepts strings of the function names. See \code{?na.action} for more details. By default it uses the na.action which is set in the global options by the respective user. The respective na.action is called on the whole data.frame supplied in the \code{data} argument. Using \code{na.omit} or \code{na.exclude} can therefore lead to unwanted results if additional unused covariates are present in the data set. By only supplying the needed data this problem can be fixed.
  }
  \item{...}{
Further arguments passed to the respective \code{surv_method}. See details.
  }
}
\details{
The primary purpose of the \code{adjustedsurv} function is to provide a convenient way to calculate confounder-adjusted survival curves using any of the methods provided in the literature. A \code{\link[=plot.adjustedsurv]{plot}} method is provided to graphically display the estimated survival curves as well. Currently the following methods can be used:

\itemize{
  \item{"\link[=surv_direct]{direct}": Direct Standardization based on a Cox-Regression model with or without stratification.}
  \item{"\link[=surv_direct_pseudo]{direct_pseudo}": Direct Standardization based on Pseudo-Values.}
  \item{"\link[=surv_iptw_km]{iptw_km}": A weighted Kaplan-Meier estimator.}
  \item{"\link[=surv_iptw_cox]{iptw_cox}": A weighted estimator based on a univariate weighted Cox-Regression model.}
  \item{"\link[=surv_iptw_pseudo]{iptw_pseudo}": A weighted estimator based on Pseudo-Values.}
  \item{"\link[=surv_matching]{matching}": Using propensity score matching to estimate the adjusted survival curves.}
  \item{"\link[=surv_emp_lik]{emp_lik}": Empirical Likelihood based estimator of the survival function.}
  \item{"\link[=surv_aiptw]{aiptw}": An Augmented Inverse Probability of Treatment Weighting estimator.}
  \item{"\link[=surv_aiptw_pseudo]{aiptw_pseudo}": An Augmented Inverse Probability of Treatment Weighting estimator using Pseudo-Values.}
  \item{"\link[=surv_tmle]{tmle}": Targeted Maximum Likelihood Estimation for survival curves.}
  \item{"\link[=surv_ostmle]{ostmle}": One-Step Targeted Maximum Likelihood Estimation for survival curves.}
  \item{"\link[=surv_km]{km}": A simple stratified Kaplan-Meier estimator without any form of adjustment.}
}

A short description of each method is contained in the documentation of the respective \code{surv_} function. For more detailed descriptions the cited literature in that same documentation can be used. The documentation for \code{method="direct"} for example can be accessed using \code{?surv_direct}.

For most methods approximations for the asymptotic variance of point estimates of the survival function have been proposed in the literature. Where available, those can be calculated and added to the output object using \code{conf_int=TRUE}. It is however recommended to use bootstrapping to estimate the variance instead, which can be done by setting \code{bootstrap=TRUE}. The \code{n_boot} argument is set to 500 by default. This number was choosen because it worked well in simulations but it does not guarantee convergence in practice. Users are recommended to inspect the bootstrapped estimates and adjust the number of replications accordingly.

To allow faster bootstrapping the user can choose to run the function on multiple CPU cores in parallel using the \code{n_cores} argument.

}
\value{
Returns an \code{adjustedsurv} object containing the following objects:

\item{adjsurv}{A \code{data.frame} of estimated adjusted survival probabilities for some points in time for each level of \code{variable}. Depending on the arguments used also includes standard errors and confidence intervals at points in time.
}
\item{method}{
The method used to adjust the survival curves.
}
\item{categorical}{Whether there are more than 2 groups in \code{variable}.
}
\item{call}{
The original function call.
}

When the argument \code{bootstrap} is set to \code{TRUE} it additionally contains the following objects:

\item{boot_data}{The adjusted survival curves estimated in each bootstrap sample.
}
\item{boot_data_same_t}{The value of the adjusted survival probabilities at the same points in time as the original survival curve for each bootstrap sample.
}
\item{boot_adjsurv}{The mean adjusted survival curves of all bootstrap samples and corresponding standard deviations and percentile confidence intervals.
}

}

\references{
Hier meine eigenen Paper
}
\author{
The function itself was written by Robin Denz, but it includes wrappers for functions written by other people. More information can be found in the respective \code{surv_method} documentation.
}

\seealso{
\code{\link{surv_aiptw}}, \code{\link{surv_aiptw_pseudo}}, \code{\link{surv_direct}}, \code{\link{surv_direct_pseudo}},
\code{\link{surv_iptw_km}}, \code{\link{surv_iptw_cox}}, \code{\link{surv_iptw_pseudo}}, \code{\link{surv_emp_lik}},
\code{\link{surv_matching}}, \code{\link{surv_tmle}}, \code{\link{surv_ostmle}}
}
\examples{
\dontrun{
library(survival)

# simulate some example data
sim_dat <- sim_confounded_surv(n=500, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# treatment assignment model
glm_mod <- glm(group ~ x2 + x3 + x5 + x6, data=sim_dat, family="binomial")

# outcome model
cox_mod <- coxph(Surv(time, event) ~ x1 + x2 + x4 + x5 + group,
                 data=sim_dat, x=TRUE)

# using direct adjustment with asymptotic confidence intervals
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="direct",
                        outcome_model=cox_mod,
                        conf_int=TRUE,
                        bootstrap=FALSE)

# using IPTW Kaplan-Meier with asymptotic confidence intervals
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="iptw_km",
                        treatment_model=glm_mod,
                        conf_int=TRUE,
                        bootstrap=FALSE)

# using AIPTW with asymptotic confidence intervals
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="aiptw",
                        outcome_model=cox_mod,
                        treatment_model=glm_mod,
                        conf_int=TRUE,
                        bootstrap=FALSE)

# using direct adjustment at custom points in time
custom_times <- c(0.001, 0.1, 0.2, 0.6, 1.1)
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="direct",
                        outcome_model=cox_mod,
                        conf_int=TRUE,
                        bootstrap=FALSE,
                        times=custom_times)

# using bootstrapping with direct adjustment
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="direct",
                        outcome_model=cox_mod,
                        conf_int=TRUE,
                        bootstrap=TRUE,
                        n_boot=500)

# using bootstrapping with direct adjustment, run in parallel
# all but one availble processor cores
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="direct",
                        outcome_model=cox_mod,
                        conf_int=TRUE,
                        bootstrap=TRUE,
                        n_boot=500,
                        n_cores=parallel::detectCores()-1)

# More specific examples can be found in the documentation of each
# respective surv_method. See ?surv_ + "method" for more information.
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
