\name{cif_direct_models}
\alias{cif_direct_models}

\title{

List of supported models in \code{cif_direct}

}
\description{

Supported models for the \code{outcome_model} argument when using \code{method="direct"} in the \code{\link[adjustedCurves]{adjustedcif}} function.

}

\details{

The following models are directly supported in the \code{outcome_model} in the \code{\link[adjustedCurves]{cif_direct}} function. The first letter in parentheses after the object name is a group indicator. Below the list there are more information for each group.

\itemize{
  \item{}{\code{\link[riskRegression]{CSC}} [\strong{A}, Required Packages: \pkg{riskRegression}]}

  \item{}{\code{\link[riskRegression]{FGR}} [\strong{B}, Required Packages: \pkg{riskRegression}]}

  \item{}{\code{\link[riskRegression]{riskRegression}} [\strong{B}, Required Packages: \pkg{riskRegression}]}

  \item{}{\code{\link[prodlim]{prodlim}} [\strong{B}, Required Packages: \pkg{prodlim}, \pkg{riskRegression}]}

  \item{}{\code{\link[randomForestSRC]{rfsrc}} [\strong{B}, Required Packages: \pkg{randomForestSRC}, \pkg{riskRegression}]}

  \item{}{\code{ARR} [\strong{B}, Required Packages: \pkg{riskRegression}]}

  \item{}{\code{\link[hal9001]{fit_hal}} [\strong{B}, Required Packages: \pkg{hal9001}, \pkg{riskRegression}]}

  \item{}{Any model with a fitting S3 prediction method or a valid \code{predict_fun} can be used as well. See below.}
}

\strong{Group A:} The direct adjusted cumulative incidences are estimated directly using the \code{\link[riskRegression]{ate}} function. Additional arguments supplied using the \code{...} syntax are passed to the \code{\link[riskRegression]{ate}} function. \cr
\strong{Group B:} The \code{\link[riskRegression]{predictRisk}} function is used to obtain predicted cumulative incidences, which are then used in the G-Computation step. Additional arguments supplied using the \code{...} syntax are passed to the \code{\link[riskRegression]{predictRisk}} function.\cr

It is sometimes possible to use models even if they are not listed here. There are two ways to make this work. The first one is to use the models S3 \code{predict} method. This works if the \code{predict} function contains the arguments \code{object}, \code{newdata}, \code{times} and \code{cause} and returns a matrix of predicted cause-specific cumulative incidences. The matrix should be of size \code{nrow(data) * length(times)}, where each row corresponds to a row in the original dataset and each column to one point in time. The matrix should contain the cause-specific cumulative incidences predicted by the model given covariates. If no such \code{predict} method exists the only option left is to write your own function which produces the output described above and supply this function to the \code{predict_fun} argument.

If you think that some important models are missing from this list, please file an issue on the official github page with a specific feature request (URL can be found in the DESCRIPTION file) or contact the package maintainer directly using the given e-mail address.

}
