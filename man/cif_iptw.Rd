\name{cif_iptw}
\alias{cif_iptw}

\title{
Inverse Probability of Treatment Weighted CIFs
}
\description{

This page explains the details of estimating inverse probability of treatment weighted cumulative incidence functions in a competing risks setting (\code{method="iptw"} in the \code{\link[=adjustedcif]{adjustedcif}} function). All regular arguments of the \code{adjustedcif} function can be used. Additionally, the \code{treatment_model} argument has to be specified in the \code{adjustedcif} call. Further arguments specific to this method are listed below.

}

\arguments{
  \item{treatment_model}{
[\strong{required}] Must be a \code{glm} or \code{multinom} model object with \code{variable} as response variable.
}
  \item{censoring_model}{
Either \code{NULL} (default) to make no adjustments for dependent censoring, or a \code{coxph} object. See \code{?ate} for more details.
  }
  \item{verbose}{
Whether to print estimation information of the \code{ate} function in the \pkg{riskRegression} package. Defaults to \code{FALSE}.
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
\item{\strong{Dependencies:} This method relies on the \pkg{riskRegression} package}
}

This method works by modeling the treatment assignment mechanism. Adjusted CIFs are calculated by first estimating appropriate case-weights for each observation in \code{data}. Those weights are used in a weighted version of the Aalen-Johansen estimator. If the weights are correctly estimated the resulting estimates will be unbiased. A more detailed description can be found in Neumann et al. (2016) and Choi et al. (2019). By utilizing another set of weights, this function can also correct the estimates for covariate-dependent censoring (Ozenne et al. 2020). Asymptotic variance calculations are based on the efficient influence curve.

Internally, this function simply calls the \code{\link[riskRegression]{ate}} function with appropriate arguments. The three-dot syntax can be used to pass further arguments to that function. It is however recommended to use the \code{\link[riskRegression]{ate}} function directly when specific settings are required.

}
\value{
Adds the following additional objects to the output of the \code{adjustedcif} function:

\itemize{
\item \code{ate_object}: The object returned by the \code{ate} function.
}
}
\references{
Anke Neumann and Cécile Billionnet (2016). "Covariate Adjustment of Cumulative Incidence Functions for Competing Risks Data Using Inverse Probability of Treatment Weighting". In: Computer Methods and Programs in Biomedicine 129, pp. 63-70

Sangbum Choi, Chaewon Kim, Hua Zhong, Eun-Seok Ryu, and Sung Won Han (2019). "Adjusted-Crude-Incidence Analysis of Multiple Treatments and Unbalanced Samples on Competing Risks". In: Statistics and Its Inference 12, pp. 423-437

Brice Maxime Hugues Ozenne, Thomas Harder Scheike, and Laila Stærk (2020). "On the Estimation of Average Treatment Effects with Right-Censored Time to Event Outcome and Competing Risks". In: Biometrical Journal 62, pp. 751-763
}
\author{
The wrapper function was written by Robin Denz, the \code{\link[riskRegression]{ate}} function itself was written by other people. See \code{?ate} for more information.
}

\seealso{
\code{\link[riskRegression]{ate}}, \code{\link{glm}}, \code{\link[nnet]{multinom}}
}
\examples{
library(adjustedCurves)

if (requireNamespace("riskRegression")) {

set.seed(42)

# simulate some data as example
sim_dat <- sim_confounded_crisk(n=50, max_t=5)
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
