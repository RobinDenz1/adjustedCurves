\name{surv_method_km}
\alias{surv_method_km}

\title{
Group-Specific Kaplan-Meier Survival Curves
}
\description{
Calculates standard Kaplan-Meier survival curves, stratified by the group variable.
No adjustment for any confounders is made. This function is included only for reference and should not be used when confounder adjusted survival curves are desired.
}
\usage{
surv_method_km(data, variable, ev_time, event, conf_int,
               conf_level=0.95, ...)
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
  \item{conf_int}{
A logical variable, indicating whether the asymptotic variances and confidence intervals of the survival probabilities should be calculated. For the available types of confidence intervals see \code{?survfit}.
  }
  \item{conf_level}{
A number specifying the confidence level of asymptotic and/or bootstrap confidence intervals.
  }
  \item{...}{
Further arguments passed to \code{survfit}.
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
The wrapper function was written by Robin Denz, the \code{survfit} function (which this wrapper is build around) was written by other people. See \code{?survfit} for more details.
}

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
