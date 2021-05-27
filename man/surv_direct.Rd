\name{surv_direct}
\alias{surv_direct}

\title{
Direct Adjusted Survival Curves
}
\description{

This page explains the details of estimating confounder-adjusted survival curves using a previously fit Cox-Regression model for single event time-to-event data using Direct Standardization. The function is called internally when using \code{method="direct"} in the \code{\link[=adjustedsurv]{adjustedsurv}} function.

}
\usage{
surv_direct(data, variable, ev_time, event, conf_int,
            conf_level=0.95, times, outcome_model, verbose=F, ...)
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
A number specifying the confidence level of asymptotic confidence intervals.
  }
  \item{outcome_model}{
Must be a \code{coxph} model object with \code{variable} as independent variable.
}
  \item{times}{
A numeric vector of time points at which the survival probability should be estimated.
}
  \item{verbose}{
Whether to print estimation information of the \code{ate} function in the \code{riskRegression} package.
}
  \item{...}{
Further arguments passed to \code{ate}.
  }
}
\details{

\itemize{
\item{\strong{Type of Adjustment:} Requires a model describing the outcome mechanism. This must be a \code{\link[survival]{coxph}} model.}
\item{\strong{Doubly-Robust:} Estimates are not Doubly-Robust.}
\item{\strong{Categorical groups:} Any number of levels in \code{variable} are allowed. Must be a factor variable.}
\item{\strong{Approximate Variance:} Calculations to approximate the variance and confidence intervals are available.}
\item{\strong{Allowed Time Values:} Allows both continuous and integer time.}
\item{\strong{Bounded Estimates:} Estimates are guaranteed to be bounded in the 0 to 1 probability range.}
\item{\strong{Monotone Function:} Estimates are guaranteed to be monotone.}
\item{\strong{Dependencies:} This method relies on the \code{\link[riskRegression]{riskRegression}} package.}
}

This method works by executing the following steps: (1) First a model is fitted which describes the outcome mechanism (time-to-event). Next (2) multiple copies of the original dataset are created, one for each possible level of the \code{variable} of interest. (3) The \code{variable} is then set to one level for all observations in each dataset. (4) The model is used to predict the survival probabilities at some points in time T for each observation in all dataset copies. (5) Those estimated probabilities are averaged for each dataset at each point in time, resulting in adjusted survival probabilities for all levels of the group variable at the specified points in time.

In the literature this method is sometimes called "Direct Standardization", "Corrected Group-Prognosis", "G-Computation" or "G-Formula". If the model in step (1) is correct this method will produce unbiased estimates of the causal survival curves. It is also one of the most efficient methods. Technically any type of model could be used. Currently only \code{\link[survival]{coxph}} models are supported however. More information can be found in the literature listed in the references.

}
\value{
Returns a \code{data.frame} object containing the estimated adjusted survival probabilities for some points in time for each level of \code{variable}.
}
\references{
I-Ming Chang, Rebecca Gelman, and Marcello Pagano. "Corrected Group Prognostic Curves and Summary Statistics". In: Journal of Chronic Diseases 35 (1982), pp. 669-674.

Robert W. Makuch. "Adjusted Survival Curve Estimation Using Covariates". In: Journal of Chronic Diseases 35.6 (1982), pp. 437-443.

Xu Zhang, Fausto R. Loberiza, John P. Klein, and Mei-Jie Zhang. "A SAS Macro for Estimation of Direct Adjusted Survival Curves Based on a Stratified Cox Regression Model". In: Computer Methods and Programs in Biomedicine 88 (2007), pp. 95-101.
}
\author{
The wrapper function was written by Robin Denz, the \code{ate} function (which this wrapper is build around) was written by other people. See \code{?ate} for more details.
}

\seealso{
\code{\link[riskRegression]{ate}}, \code{\link[survival]{coxph}}
}
\examples{
\dontrun{
library(survival)

# simulate some data as example
sim_dat <- sim_confounded_surv(n=500, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# estimate a cox-regression for the outcome
cox_mod <- coxph(Surv(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6 + group,
                 data=sim_dat, x=TRUE)


# use it to calculate adjusted survival curves
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="direct",
                        outcome_model=cox_mod,
                        conf_int=FALSE)

# plot the curves
plot(adjsurv)
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
