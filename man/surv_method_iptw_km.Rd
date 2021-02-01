\name{surv_method_iptw_km}
\alias{surv_method_iptw_km}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Inverse Probability of Treatment Weighted Kaplan-Meier estimates
}
\description{

}
\usage{
surv_method_iptw_km(data, variable, ev_time, event, sd,
                    treatment_model, weight_method="ps", ...)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{weight_method}{
Method used in \code{WeightIt} function call. Ignored if \code{treatment_model} is not a formula object.
}
  \item{...}{
Further arguments passed to \code{WeightIt}.
  }
}
\details{
BLA
}
\value{
Returns a \code{data.frame} object containing the estimated adjusted survival probabilities for some points in time for each level of \code{variable}.
}
\references{
Jun Xie and Chaofeng Liu. "Adjusted Kaplan-Meier Estimator and Log-Rank Test with Inverse Probability of Treatment Weighting for Survival Data". In: Statistics in Medicine 24 (2005), pp. 3089-3110.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{WeightIt}}, \code{RISCA}
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
