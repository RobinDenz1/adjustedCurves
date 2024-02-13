\name{models_surv_direct}
\alias{models_surv_direct}

\title{

List of supported models in \code{surv_direct}

}
\description{

Supported models for the \code{outcome_model} argument when using \code{method="direct"} in the \code{\link[adjustedCurves]{adjustedsurv}} function.

}

\details{

The following models are directly supported in the \code{outcome_model} in the \code{\link[adjustedCurves]{surv_direct}} function. The first letter in parentheses after the object name is a group indicator. Below the list there are more information for each group.

\itemize{
  \item \code{\link[survival]{coxph}} [\strong{A}, Required Packages: \pkg{survival}, \pkg{riskRegression}]

  \item \code{\link[rms]{cph}} [\strong{A}, Required Packages: \pkg{rms}, \pkg{survival}, \pkg{riskRegression}]

  \item \code{\link[timereg]{aalen}} [\strong{B}, Required Packages: \pkg{timereg}, \pkg{pec}]

  \item \code{\link[timereg]{cox.aalen}} [\strong{B}, Required Packages: \pkg{timereg}, \pkg{pec}]

  \item \code{\link[riskRegression]{selectCox}} [\strong{B}, Required Packages: \pkg{riskRegression}, \pkg{pec}]

  \item \code{\link[pec]{pecCforest}} [\strong{B}, Required Packages: \pkg{pec}]

  \item \code{\link[pec]{pecRpart}} [\strong{B}, Required Packages: \pkg{pec}, Bootstrapping not allowed.]

  \item \code{\link[riskRegression]{riskRegression}} [\strong{C}, Required Packages: \pkg{riskRegression}]

  \item \code{\link[prodlim]{prodlim}} [\strong{C}, Required Packages: \pkg{prodlim}, \pkg{riskRegression}]

  \item \code{\link[rms]{psm}} [\strong{C}, Required Packages: \pkg{rms}, \pkg{riskRegression}]

  \item \code{flexsurvreg} [\strong{C}, Required Packages: \pkg{flexsurv}, \pkg{riskRegression}]

  \item \code{flexsurvspline} [\strong{C}, Required Packages: \pkg{flexsurv}, \pkg{riskRegression}]

  \item \code{\link[ranger]{ranger}} [\strong{C}, Required Packages: \pkg{ranger}, \pkg{riskRegression}]

  \item \code{rfsrc} [\strong{C}, Required Packages: \pkg{randomForestSRC}, \pkg{riskRegression}]

  \item \code{ARR} [\strong{C}, Required Packages: \pkg{riskRegression}]

  \item \code{\link[riskRegression]{penalizedS3}} [\strong{C}, Required Packages: \pkg{penalized}, \pkg{riskRegression}]

  \item \code{gbm} [\strong{C}, Required Packages: \pkg{gbm}, \pkg{riskRegression}]

  \item \code{fit_hal} [\strong{C}, Required Packages: \pkg{hal9001}, \pkg{riskRegression}]

  \item \code{fitSmoothHazard} [\strong{C}, Required Packages: \pkg{casebase}, \pkg{riskRegression}]

  \item \code{\link[stats]{glm}} [\strong{D}, Required Packages: \pkg{stats}, \pkg{pec}]

  \item \code{\link[rms]{ols}} [\strong{D}, Required Packages: \pkg{rms}, \pkg{pec}]

  \item \code{randomForest} [\strong{D}, Required Packages: \pkg{randomForest}, \pkg{pec}]

  \item \code{mexhaz} [\strong{E}, Required Packages: \pkg{mexhaz}]

  \item Any model with a fitting S3 prediction method or a valid \code{predict_fun} can be used as well. See below.
}

\strong{Group A:} The direct adjusted survival probabilities are estimated directly using the \code{\link[riskRegression]{ate}} function. Additional arguments supplied using the \code{...} syntax are passed to the \code{\link[riskRegression]{ate}} function. Note that \code{Surv()} calls required to fit the model should be made inside the formula, not saved elsewhere. \cr
\strong{Group B:} Predicted survival probabilities are obtained using the \code{\link[pec]{predictSurvProb}} function. The G-Computation is carried out using those. Additional arguments supplied using the \code{...} syntax are passed to the \code{\link[pec]{predictSurvProb}} function. \cr
\strong{Group C:} The \code{\link[riskRegression]{predictRisk}} function is used to obtain predicted cumulative incidences, which are then transformed to survival probabilities. Additional arguments supplied using the \code{...} syntax are passed to the \code{\link[riskRegression]{predictRisk}} function.\cr
\strong{Group D:} These models are only allowed if there is no censoring. Predicted survival probabilities are obtained using the \code{predictProb} function from the \pkg{pec} package. Additional arguments supplied using the \code{...} syntax are passed to the \code{predictProb} function.\cr
\strong{Group E:} Custom code is used to obtain predicted survival probabilities. Additional arguments are not used.

It is sometimes possible to use models even if they are not listed here. There are two ways to make this work. The first one is to use the models S3 \code{predict} method. This works if the \code{predict} function contains the arguments \code{object}, \code{newdata} and \code{times} and returns a matrix of predicted survival probabilities. The matrix should be of size \code{nrow(data) * length(times)}, where each row corresponds to a row in the original dataset and each column to one point in time. The matrix should contain the survival probabilities predicted by the model given covariates. If no such \code{predict} method exists the only option left is to write your own function which produces the output described above and supply this function to the \code{predict_fun} argument.

If you think that some important models are missing from this list, please file an issue on the official github page with a specific feature request (URL can be found in the DESCRIPTION file) or contact the package maintainer directly using the given e-mail address.

}

\note{
When using outcome models which are not directly supported (either through the default predict method or a custom \code{predict_fun}) it might be necessary to set the \code{clean_data} argument of the \code{adjustedsurv} function to \code{FALSE}.
}
