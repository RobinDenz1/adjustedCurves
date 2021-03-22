\name{cif_aiptw_pseudo}
\alias{cif_aiptw_pseudo}

\title{
Augmented Inverse Probability of Treatment Weighted CIFs using Pseudo-Values
}
\description{
Function to calculate augmented inverse probability of treatment weighted CIFs using pseudo-values in a competing risks setting. Instead of only modeling the outcome mechanism or the treatment assignment mechanism, both kind of models are required to use this method. If either of those models are correctly specified, unbiased estimates will be obtained. Both the outcome and treatment assignment part are modeled using pseudo-values. See details for more information.
}
\usage{
cif_aiptw_pseudo(data, variable, ev_time, event, cause, conf_int,
                 conf_level=0.95, times, outcome_vars,
                 treatment_model, type_time="factor", spline_df=10)
}

\arguments{
  \item{data}{
A \code{data.frame} object containing the needed time-to-event data in standard format.
}
  \item{variable}{
A character string specifying the variable by which the CIFs should be grouped. Must be a valid columname of \code{data}.
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
  \item{conf_int}{
A logical variable, indicating whether the asymptotic variances and confidence intervals of the CIFs should be calculated.
  }
  \item{conf_level}{
A number specifying the confidence level of asymptotic and/or bootstrap confidence intervals.
  }
  \item{times}{
A numeric vector of time points at which the CIFs should be estimated.
}
  \item{outcome_vars}{
A character vector of column names specifying variables to be used when modeling the outcome mechanism using \code{geese}. See details and examples.
}
  \item{treatment_model}{
Must be a \code{glm} or \code{multinom} model object with \code{variable} as response variable. See details and examples.
}
  \item{type_time}{
A character string specifying how the time should be modeled. See details.
}
  \item{spline_df}{
The number of degrees of freedom used for the natural-spline or B-spline function. Ignored if \code{type_time="factor"}.
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
\item{\strong{Dependencies:} This method relies on the \code{geepack} and \code{\link[prodlim]{prodlim}} packages.}
}

Instead of only modeling the outcome mechanism or the treatment assignment mechanism, both kind of models are required to use this method. If either of those models are correctly specified, unbiased estimates will be obtained. In contrast to the \code{\link[=cif_aiptw]{"aiptw"}} method, the "aiptw_pseudo" method uses a generalised estimation equation (geese) approach to model the outcome mechanism. The model is fit in the same way as described in the \code{\link[=cif_direct_pseudo]{"direct_pseudo"}} method. Those Direct Standardization based estimates are then transformed using the previously estimated propensity score. This results in the doubly-robust property of the method. More information on this particular method can be found in the original article by Wang (2018). The original article only deals with survival probabilities without competing risks, but the only difference to the CIF estimation with competing risks is the calculation of the pseudo-values. More information on Pseudo-Values is available in Andersen et al. (2017) and Andersen and Perme (2010).

When estimating the \code{geese} model the \code{ev_time} variable is used as a factor by default. This results in one coefficient being estimated for each unique point in time, which can be very slow computationally if there are a lot of unique points in time and/or the dataset has many rows. In these cases it is recommended to use \code{type_time="bs"} or \code{type_time="ns"}, which results in the \code{ev_time} being modeled using B-Splines or Natural Splines. Simulation studies indicate that there is little difference in the estimates when an appropriately large number of \code{spline_df} is used.

}
\value{
Returns a \code{data.frame} object containing the estimated adjusted survival probabilities for some points in time for each level of \code{variable}.
}
\references{

Jixian Wang. "A Simple, Doubly Robust, Efficient Estimator for Survival Functions Using Pseudo Observations". In: Pharmaceutical Statistics 17.38-48 (2018).

James M. Robins and Andrea Rotnitzky. "Recovery of Information and Adjustment for Dependent Censoring Using Surrogate Markers". In: AIDS Epidemiology: Methodological Issues. Ed. by Nicholas P. Jewell, Klaus Dietz, and Vernon T. Farewell. New York: Springer Science + Business Media, 1992, pp. 297-331.

Per K. Andersen, Elisavet Syriopoulou, and Erik T. Parner. "Causal Inference in Survival Analysis using Pseudo-Observations". In: Statistics in Medicine 36 (2017), pp. 2669-2681.

Per K. Andersen and Maja Pohar Perme. "Pseudo-Observations in Survival Analysis". In: Statistical Methods in Medical Research 19 (2010), pp. 71-99.

Aris Perperoglou, Willi Sauerbrei, Michal Abrahamowicz, and Matthias Schmid. "A Review of Spline Function Procedures in R". In: BMC Medical Research Methodology 19 (2019), pp. 1-16
}
\author{
Robin Denz
}

\seealso{
\code{\link[geepack]{geese}}, \code{\link[prodlim]{jackknife}}, \code{\link[splines]{ns}}, \code{\link[splines]{bs}}
}
\examples{
\dontrun{
# simulate some data as example
sim_dat <- sim_confounded_surv(n=500, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# estimate a treatment assignment model
glm_mod <- glm(group ~ x1 + x3 + x5 + x6, data=sim_dat, family="binomial")

# use it + pseudo values + geese model to calculate adjusted survival curves
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="aiptw_pseudo",
                        outcome_vars=c("x1", "x2", "x3", "x4", "x5", "x6"),
                        treatment_model=glm_mod,
                        conf_int=TRUE)

# plot the curves
plot(adjsurv, draw_ci=TRUE)
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
