\name{surv_method_aiptw}
\alias{surv_method_aiptw}

\title{
Augmented Inverse Probability of Treatment Weighted Survival Curves
}
\description{
Function to calculate augmented inverse probability of treatment weighted survival curves. Instead of only modeling the outcome mechanism or the treatment assignment mechanism, both kind of models are required to use this method. If either of those models are correctly specified, unbiased estimates will be obtained. Can also be used to adjust for dependent censoring. See details for more information.
}
\usage{
surv_method_aiptw(data, variable, ev_time, event, sd, times,
                  outcome_model=NULL, treatment_model=NULL,
                  censoring_model=NULL, verbose=F, ...)
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
  \item{sd}{
A logical variable, indicating whether the standard deviation of the survival probabilities should be calculated. Needed to calculate confidence intervals.
}
  \item{times}{
A numeric vector of time points at which the survival probability should be estimated or \code{NULL}. If \code{NULL} the survival is estimated at all points in time at which an event occured.
}
  \item{outcome_model}{
Must be a \code{coxph} model object, modeling the time-to-event mechanism. See details and examples.
}
  \item{treatment_model}{
Must be a \code{glm} model object with \code{variable} as response variable. See details and examples.
}
  \item{censoring_model}{
Must be a \code{coxph} model object, modeling the censoring mechanism. See details and examples.
}
  \item{verbose}{
Whether to print estimation information of the \code{ate} function in the \code{riskRegression} package.
}
  \item{...}{
Further arguments passed to \code{ate}.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a \code{data.frame} object containing the estimated adjusted survival probabilities for some points in time for each level of \code{variable}.
}
\references{

James M. Robins and Andrea Rotnitzky. "Recovery of Information and Adjustment for Dependent Censoring Using Surrogate Markers". In: AIDS Epidemiology: Methodological Issues. Ed. by Nicholas P. Jewell, Klaus Dietz, and Vernon T. Farewell. New York: Springer Science + Business Media, 1992, pp. 297-331.

Alan E. Hubbard, Mark J. van der Laan, and James M. Robins. "Nonparametric Locally Efficient Estimation of the Treatment Specific Survival Distribution with Right Censored Data and Covariates in Observational Studies". In: Statistical Models in Epidemiology, the Environment, and Clinical Trials. Ed. by M. Elizabeth Halloran and Donald Berry. New York: Springer Science + Business Media, 2000, pp. 135-177.

Min Zhang and Douglas E. Schaubel. "Contrasting Treatment-Specific Survival Using Double-Robust Estimators". In: Statistics in Medicine 31.30 (2012), pp. 4255-4268.

Xiaofei Bai, Anastasios A. Tsiatis, and Sean M. O'Brien. "Doubly-Robust Estimators of Treatment-Specific Survival Distributions in Observational Studies with Stratified Sampling". In: Biometrics 69 (2013), pp. 830-839.

Brice Maxime Hugues Ozenne, Thomas Harder Scheike, and Laila Staerk. "On the Estimation of Average Treatment Effects with Right-Censored Time to Event Outcome and Competing Risks". In: Biometrical Journal 62 (2020), pp. 751-763.

}
\author{
The wrapper function was written by Robin Denz, the \code{ate} function (which this wrapper is build around) was written by other people. See \code{?ate} for more details.
}

\seealso{
\code{\link{ate}}, \code{\link{coxph}}, \code{\link{glm}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
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
