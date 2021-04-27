\name{cif_iptw}
\alias{cif_iptw}

\title{
Inverse Probability of Treatment Weighted CIFs
}
\description{
This method works by modelling the treatment assignment mechanism. The individual observations are weighted using inverse probability of treatment weights or other balancing weights (see details). If the model is correctly specified, unbiased results will be obtained. See details for more information.
}
\usage{
cif_iptw(data, variable, ev_time, event, cause, conf_int,
         conf_level=0.95, times, treatment_model,
         censoring_model=NULL, verbose=F, ...)
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
A logical variable, indicating whether the asymptotic variances and confidence intervals of the CIFs should be calculated. Not available for all methods.
  }
  \item{conf_level}{
A number specifying the confidence level of asymptotic confidence intervals.
  }
  \item{times}{
A numeric vector of time points at which the cumulative incidences should be estimated.
  }
  \item{treatment_model}{
Must be a \code{glm} or \code{multinom} model object with \code{variable} as response variable.
}
  \item{censoring_model}{
Either \code{NULL} to make no adjustments for dependent censoring, or a \code{coxph} object. See \code{?ate} for more details.
  }
  \item{verbose}{
Whether to print estimation information of the \code{ate} function in the \code{riskRegression} package.
  }
  \item{...}{
Further arguments passed to \code{\link[riskRegression]{ate}}.
  }
}
\details{

\itemize{
\item{\strong{Type of Adjustment:} Requires a model describing the treatment assignment mechanism. This must be either a \code{\link{glm}} or a \code{\link{multinom}} object.}
\item{\strong{Doubly-Robust:} Estimates are not Doubly-Robust.}
\item{\strong{Categorical groups:} Any number of levels in \code{variable} are allowed. Must be a factor variable.}
\item{\strong{Approximate Variance:} Calculations to approximate the variance and confidence intervals are available.}
\item{\strong{Allowed Time Values:} Allows both continuous and integer time.}
\item{\strong{Bounded Estimates:} Estimates are guaranteed to be bounded in the 0 to 1 probability range.}
\item{\strong{Monotone Function:} Estimates are guaranteed to be monotone.}
\item{\strong{Dependencies:} This method relies on the \code{\link[riskRegression]{riskRegression}} package}
}

Adjusted CIFs are calculated by first estimating appropriate case-weights for each observation in \code{data}. Those weights are used in a weighted version of the Aalen-Johansen estimator. If the weights are correctly estimated the resulting estimates will be unbiased. A more detailed description can be found in Neumann et al. (2016) and Choi et al. (2019). By utilizing another set of weights, this function can also correct the estimates for covariate-dependent censoring (Ozenne et al. 2020). Asymptotic variance calculations are based on the efficient influence curve.

Internally, this function simply calls the \code{\link[riskRegression]{ate}} function with appropriate arguments. The three-dot syntax can be used to pass further arguments to that function. It is however recommened to use the \code{\link[riskRegression]{ate}} function directly when specific settings are required.

}
\value{
Returns a \code{data.frame} object containing the estimated adjusted CIFs for some points in time for each level of \code{variable}.
}
\references{
Anke Neumann and Cecile Billionnet. "Covariate Adjustment of Cumulative Incidence Functions for Competing Risks Data Using Inverse Probability of Treatment Weighting". In: Computer Methods and Programs in Biomedicine 129 (2016), pp. 63–70.

Sangbum Choi, Chaewon Kim, Hua Zhong, Eun-Seok Ryu, and SungWon Han. "Adjusted-Crude-Incidence Analysis of Multiple Treatments and Unbalanced Samples on Competing Risks". In: Statistics and Its Inference 12 (2019), pp. 423–437.

Brice Maxime Hugues Ozenne, Thomas Harder Scheike, and Laila Staerk. "On the Estimation of Average Treatment Effects with Right-Censored Time to Event Outcome and Competing Risks". In: Biometrical Journal 62 (2020), pp. 751–763.
}
\author{
The wrapper function was written by Robin Denz, the \code{\link[riskRegression]{ate}} function itself was written by other people. See \code{?ate} for more information.
}

\seealso{
\code{\link[riskRegression]{ate}}, \code{\link{glm}}, \code{\link[nnet]{multinom}}
}
\examples{
\dontrun{
# simulate some data as example
sim_dat <- sim_confounded_crisk(n=500, max_t=5)
sim_dat$group <- as.factor(sim_dat$group)

# estimate a treatment assignment model
glm_mod <- glm(group ~ x1 + x3 + x5 + x6, data=sim_dat, family="binomial")

# use it to calculate adjusted CIFs
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      cause=1,
                      method="iptw",
                      treatment_model=glm_mod)
plot(adjcif)
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
