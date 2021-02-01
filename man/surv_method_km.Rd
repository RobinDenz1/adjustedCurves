\name{surv_method_km}
\alias{surv_method_km}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Group-Specific Kaplan-Meier Survival Curves
}
\description{
Calculates standard Kaplan-Meier survival curves, stratified by the group variable.
No adjustment for any confounders is made. This function is included only for reference and should not be used when confounder adjusted survival curves are desired.
}
\usage{
surv_method_km(data, variable, ev_time, event, sd, ...)
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
  \item{...}{
Further arguments passed to \code{survfit()}.
  }
}
\details{

}
\value{
Returns a \code{data.frame} object containing the estimated survival probabilities for some points in time for each level of \code{variable}.
}
\references{
E. L. Kaplan and Paul Meier. "Nonparametric Estimation from Incomplete Observations". In: Journal of the American Statistical Association 53.282 (1958), pp. 457-481.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{survival}}, \code{\link{survfit}}
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
