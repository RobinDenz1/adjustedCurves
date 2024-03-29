\name{surv_direct}
\alias{surv_direct}

\title{
Direct Adjusted Survival Curves
}
\description{

This page explains the details of estimating confounder-adjusted survival curves using a previously fit Cox-Regression model for single event time-to-event data using Direct Standardization (\code{method="direct"} in the \code{\link[=adjustedsurv]{adjustedsurv}} function). All regular arguments of the \code{adjustedsurv} function can be used. Additionally, the \code{outcome_model} argument has to be specified in the \code{adjustedsurv} call. Further arguments specific to this method are listed below.

}
\arguments{
  \item{outcome_model}{
[\strong{required}] Must be a previously fit model object including \code{variable} as independent variable. Apart from the classic \code{coxph} model this function also supports a variety of other models. See \code{\link[adjustedCurves]{models_surv_direct}} for a list of supported model objects and some more details.
}
  \item{verbose}{
Whether to print estimation information of the \code{ate} function in the \pkg{riskRegression} package. Ignored if \code{outcome_model} is not a \code{coxph} object. Defaults to \code{FALSE}.
}
  \item{predict_fun}{
A function which should be used to calculate the predicted survival probabilities given covariates and some points in time. This argument only needs to be specified if the kind of model supplied in the \code{outcome_model} is not directly supported. See \code{\link[adjustedCurves]{models_surv_direct}} for more information. Defaults to \code{NULL}.
}
  \item{...}{
Further arguments passed to \code{ate} if \code{outcome_model} is a \code{coxph} object. Otherwise the additional arguments are passed to the respective \code{predict} method. See \code{\link[adjustedCurves]{models_surv_direct}} for more information.
  }
}
\details{

\itemize{
\item{\strong{Type of Adjustment:} Requires a model describing the outcome mechanism. See \code{\link[adjustedCurves]{models_surv_direct}} for a list of supported model objects and some more details.}
\item{\strong{Doubly-Robust:} Estimates are not Doubly-Robust.}
\item{\strong{Categorical groups:} Any number of levels in \code{variable} are allowed. Must be a factor variable.}
\item{\strong{Approximate Variance:} Calculations to approximate the variance and confidence intervals are available only if \code{outcome_model} is a \code{coxph} object. The \code{\link[riskRegression]{ate}} function is used for the calculation in that case. Bootstrap confidence intervals can however be calculated with all supported models. See \code{?adjustedsurv} for more information on bootstrapping.}
\item{\strong{Allowed Time Values:} Allows both continuous and integer time.}
\item{\strong{Bounded Estimates:} Estimates are guaranteed to be bounded in the 0 to 1 probability range.}
\item{\strong{Monotone Function:} Estimates are guaranteed to be monotone.}
\item{\strong{Dependencies:} This method relies on the \pkg{riskRegression} package. Depending on \code{outcome_model} other packages might be needed. See \code{\link[adjustedCurves]{models_surv_direct}} for more details.}
}

This method works by executing the following steps: (1) First a model is fitted which describes the outcome mechanism (time-to-event). Next (2) multiple copies of the original dataset are created, one for each possible level of the \code{variable} of interest. (3) The \code{variable} is then set to one level for all observations in each dataset. (4) The model is used to predict the survival probabilities at some points in time T for each observation in all dataset copies. (5) Those estimated probabilities are averaged for each dataset at each point in time, resulting in adjusted survival probabilities for all levels of the group variable at the specified points in time.

In the literature this method is sometimes called "Direct Standardization", "Corrected Group-Prognosis", "G-Computation" or "G-Formula". If the model in step (1) is "correct"" this method will produce unbiased estimates of the counterfactual survival curves. A model can be called a "correct" model in this context if it can be used to produce unbiased estimates of the true (but unknown) individual survival probabilities given covariates. When used properly this is one of the most efficient methods. More information can be found in the literature listed in the references. The most popular model for describing the outcome mechanism in a time-to-event context is the Cox-regression model (\code{\link[survival]{coxph}}). This function however also supports a variety of other models.

}
\value{
Adds the following additional objects to the output of the \code{adjustedsurv} function:

\itemize{
\item \code{ate_object}: The object returned by the \code{ate} function.
}

}
\references{
I-Ming Chang, Rebecca Gelman, and Marcello Pagano (1982). "Corrected Group Prognostic Curves and Summary Statistics". In: Journal of Chronic Diseases 35, pp. 669-674

Robert W. Makuch (1982). "Adjusted Survival Curve Estimation Using Covariates". In: Journal of Chronic Diseases 35.6, pp. 437-443

Xu Zhang, Fausto R. Loberiza, John P. Klein, and Mei-Jie Zhang (2007). "A SAS Macro for Estimation of Direct Adjusted Survival Curves Based on a Stratified Cox Regression Model". In: Computer Methods and Programs in Biomedicine 88, pp. 95-101
}
\author{
The function itself was written by Robin Denz. When using \code{coxph} models however, this function is just a wrapper around the \code{\link[riskRegression]{ate}} function, which was written by other people. See \code{?ate} for more information.
}

\seealso{
\code{\link[adjustedCurves]{models_surv_direct}}, \code{\link[riskRegression]{ate}}, \code{\link[survival]{coxph}}
}
\examples{
library(adjustedCurves)
library(survival)

if (requireNamespace("riskRegression")) {

library(riskRegression)

set.seed(42)

# simulate some data as example
sim_dat <- sim_confounded_surv(n=50, max_t=1.2)
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

# not run to avoid dependency on flexsurv and mice too slow
if (interactive()) {
## using a flexsurv() model, this requires the 'fleysurv' package
mod_flexsurvreg <- flexsurvreg(Surv(time, event) ~ group + x1 + x2 + x5 + x6,
                               data=sim_dat, dist="gengamma")

# using it to calculate the adjusted survival curves
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="direct",
                        outcome_model=mod_flexsurvreg,
                        conf_int=FALSE)

# plot using steps=FALSE to draw them as smooth functions, since
# they were estimated using a parametric model
plot(adjsurv, steps=FALSE)
}

\donttest{
## using multiple imputation
if (requireNamespace("mice")) {
library(mice)

# introduce random missingness in x1 as example
# NOTE: This is only done as an example, in reality you would
#       already have missing data, not introduce it yourself.
sim_dat$x1 <- ifelse(runif(n=50) < 0.5, sim_dat$x1, NA)

# perform multiple imputation
mids <- mice::mice(data=sim_dat, method="pmm", m=5, printFlag=FALSE)

# fit model for each imputed dataset
mira <- with(mids, coxph(Surv(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6 + group,
                         x=TRUE))

# calculate adjusted survival curves on imputed data
adj <- adjustedsurv(data=mids,
                    variable="group",
                    ev_time="time",
                    event="event",
                    method="direct",
                    outcome_model=mira)
plot(adj)
}
}
}
}
