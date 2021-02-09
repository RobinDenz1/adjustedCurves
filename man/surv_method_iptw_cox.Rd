\name{surv_method_iptw_cox}
\alias{surv_method_iptw_cox}

\title{
Inverse Probability of Treatment Weighted Survival using Cox-Regression
}
\description{
This method works by modelling the treatment assignment mechanism. The individual observations are weighted using inverse probability of treatment weights or other balancing weights (see details) and a univariate cox-regression model. If the model is correctly specified, unbiased results will be obtained. The only difference to the \code{surv_method_iptw_km} is a slightly different weighting approach. They are essentially equal in large samples. See details for more information.
}
\usage{
surv_method_iptw_cox(data, variable, ev_time, event, sd,
                     treatment_model, weight_method="ps", ...)
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
  \item{weight_method}{
Method used in \code{WeightIt} function call. Ignored if \code{treatment_model} is not a formula object.
}
  \item{...}{
Further arguments passed to \code{WeightIt}.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a \code{data.frame} object containing the estimated adjusted survival probabilities for some points in time for each level of \code{variable}.
}
\references{
Stephen R. Cole and Miguel A. Hernan. "Adjusted Survival Curves with Inverse Probability Weights". In: Computer Methods and Programs in Biomedicine 2003.75 (2004), pp. 45-49.
}
\author{
Robin Denz
}

\seealso{
\code{\link{WeightIt}}, \code{\link{coxph}}
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
