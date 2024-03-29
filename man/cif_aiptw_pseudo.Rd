\name{cif_aiptw_pseudo}
\alias{cif_aiptw_pseudo}

\title{
Augmented Inverse Probability of Treatment Weighted CIFs using Pseudo-Values
}
\description{

This page explains the details of estimating augmented inverse probability of treatment weighted CIFs using pseudo-values in a competing risks setting (\code{method="aiptw_pseudo"} in the \code{\link[=adjustedcif]{adjustedcif}} function). All regular arguments of the \code{adjustedcif} function can be used. Additionally, the \code{outcome_vars} argument and the \code{treatment_model} argument have to be specified in the \code{adjustedcif} call. Further arguments specific to this method are listed below.

}

\arguments{
  \item{outcome_vars}{
[\strong{required}] A character vector of column names specifying variables to be used when modeling the outcome mechanism using \code{geese}. See details and examples.
}
  \item{treatment_model}{
[\strong{required}] Must be a \code{glm} or \code{multinom} model object with \code{variable} as response variable. See details and examples.
}
  \item{type_time}{
A character string specifying how the time should be modeled. Possible values are \code{"factor"} (modeling each point in time as a separate variable, the default), \code{"bs"} (modeling time using B-Splines) or \code{"ns"} (modeling time using natural splines).
}
  \item{spline_df}{
The number of degrees of freedom used for the natural-spline or B-spline function. Defaults to 5. Ignored if \code{type_time="factor"}.
}
}
\details{

\itemize{
\item{\strong{Type of Adjustment:} Requires a treatment assignment model (\code{\link{glm}} or \code{\link{multinom}}) and a character vector of variable names used to model the outcome mechanism (internally uses \code{\link[geepack]{geese}}). In contrast to the \code{\link[=cif_aiptw]{"aiptw"}} method this function does not allow for dependent censoring.}
\item{\strong{Doubly-Robust:} Estimates are Doubly-Robust.}
\item{\strong{Categorical groups:} Any number of levels in \code{variable} are allowed. Must be a factor variable.}
\item{\strong{Approximate Variance:} Calculations to approximate the variance and confidence intervals are available.}
\item{\strong{Allowed Time Values:} Allows both continuous and integer time.}
\item{\strong{Bounded Estimates:} Estimates are not guaranteed to be bounded in the 0 to 1 probability range.}
\item{\strong{Monotone Function:} Estimates are not guaranteed to be monotone.}
\item{\strong{Dependencies:} This method relies on the \pkg{geepack} and \pkg{prodlim} packages.}
}

Instead of only modeling the outcome mechanism or the treatment assignment mechanism, both kind of models are required to use this method. If either of those models are correctly specified, unbiased estimates will be obtained. In contrast to the \code{\link[=cif_aiptw]{"aiptw"}} method, the "aiptw_pseudo" method uses a generalized estimation equation (geese) approach to model the outcome mechanism. The model is fit in the same way as described in the \code{\link[=cif_direct_pseudo]{"direct_pseudo"}} method. Those Direct Standardization based estimates are then transformed using the previously estimated propensity score. This results in the doubly-robust property of the method. More information on this particular method can be found in the original article by Wang (2018). The original article only deals with survival probabilities without competing risks, but the only difference to the CIF estimation with competing risks is the calculation of the pseudo-values. More information on Pseudo-Values is available in Andersen et al. (2017) and Andersen and Perme (2010).

When estimating the \code{geese} model the \code{ev_time} variable is used as a factor by default. This results in one coefficient being estimated for each unique point in time, which can be very slow computationally if there are a lot of unique points in time and/or the dataset has many rows. In these cases it is recommended to use \code{type_time="bs"} or \code{type_time="ns"}, which results in the \code{ev_time} being modeled using B-Splines or Natural Splines. Simulation studies indicate that there is little difference in the estimates when an appropriately large number of \code{spline_df} is used.

}
\value{
Adds the following additional objects to the output of the \code{adjustedcif} function:

\itemize{
\item \code{pseudo_values}: The matrix of estimated pseudo-values.
\item \code{geese_model}: The geese model used to make the predictions.
}
}
\references{
Jixian Wang (2018). "A Simple, Doubly Robust, Efficient Estimator for Survival Functions Using Pseudo Observations". In: Pharmaceutical Statistics 17.38-48

James M. Robins and Andrea Rotnitzky (1992). "Recovery of Information and Adjustment for Dependent Censoring Using Surrogate Markers". In: AIDS Epidemiology: Methodological Issues. Ed. by Nicholas P. Jewell, Klaus Dietz, and Vernon T. Farewell. New York: Springer Science + Business Media, pp. 297-331

Per Kragh Andersen, Elisavet Syriopoulou, and Erik T. Parner (2017). "Causal Inference in Survival Analysis using Pseudo-Observations". In: Statistics in Medicine 36, pp. 2669-2681

Per Kragh Andersen and Maja Pohar Perme (2010). "Pseudo-Observations in Survival Analysis". In: Statistical Methods in Medical Research 19, pp. 71-99

Aris Perperoglou, Willi Sauerbrei, Michal Abrahamowicz, and Matthias Schmid (2019). "A Review of Spline Function Procedures in R". in: BMC Medical Research Methodology 19.46, pp. 1-16
}
\author{

Jixian Wang supplied the R source code used in the original article, which was used by Robin Denz to create a generalized version of this method with additional functionality and improved performance.

}

\seealso{
\code{\link[geepack]{geese}}, \code{\link[prodlim]{jackknife}}, \code{\link[splines]{ns}}, \code{\link[splines]{bs}}
}
\examples{
library(adjustedCurves)

if (requireNamespace("geepack") & requireNamespace("prodlim")) {

set.seed(42)

# simulate some data as example
sim_dat <- sim_confounded_crisk(n=30, max_t=5)
sim_dat$group <- as.factor(sim_dat$group)

# estimate a treatment assignment model
glm_mod <- glm(group ~ x1 + x3 + x5 + x6, data=sim_dat, family="binomial")

# use it + pseudo values + geese model to calculate adjusted CIFs
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      cause=1,
                      method="aiptw_pseudo",
                      outcome_vars=c("x1", "x2", "x3", "x4", "x5", "x6"),
                      treatment_model=glm_mod,
                      conf_int=FALSE,
                      force_bounds=TRUE,
                      iso_reg=TRUE)

# plot the curves
plot(adjcif)
}
}
