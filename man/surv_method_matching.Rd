\name{surv_method_matching}
\alias{surv_method_matching}

\title{
Using Matching to Calculate Adjusted Survival Curves
}
\description{
This method works by first matching controls and treated patients using a previously estimated propensity score. Using this new matched population, a simple Kaplan-Meier or Nelson-Aalen estimator is used to calculate the adjusted survival curves. Simulation studies have shown that this method, while unbiased if used appropriately, is significantly less efficient than other methods. See details.
}
\usage{
surv_method_matching(data, variable, ev_time, event, sd,
                     treatment_model, ...)
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
  \item{treatment_model}{
Must be either a model object with \code{variable} as response variable, a vector of weights or a formula which can be passed to \code{WeightIt}.
}
  \item{...}{
Further arguments passed to the \code{Match} function of the \code{Matching} Package.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a \code{data.frame} object containing the estimated adjusted survival probabilities for some points in time for each level of \code{variable}.
}
\references{
Angela Winnett and Peter Sasieni. "Adjusted Nelson-Aalen Estimates with Retrospective Matching". In: Journal of the American Statistical Association 97.457 (2002), pp. 245-256.

Peter C. Austin, Neal Thomas, and Donald B. Rubin. "Covariate-Adjusted Survival Analyses in Propensity-Score Matched Samples: Imputing Potential Time-To-Event Outcomes". In: Statistical Methods in Medical Research 29.3 (2020), pp. 728-751.
}
\author{
The wrapper function was written by Robin Denz, the \code{Match} function (which this wrapper is build around) was written by other people. See \code{?Match} for more details.
}

\seealso{
\code{\link[Matching]{Match}}, \code{\link[survival]{survfit}}
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
