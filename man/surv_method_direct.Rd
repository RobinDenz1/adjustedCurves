\name{surv_method_direct}
\alias{surv_method_direct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Direct Adjusted Survival Curves using Cox-Regression
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
surv_method_direct(data, variable, ev_time, event, sd,
                   outcome_model, times, verbose=F, ...)
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
  \item{outcome_model}{
Must be either a model object with \code{variable} as response variable, a vector of weights or a formula which can be passed to \code{WeightIt}.
}
  \item{times}{
A numeric vector of time points at which the survival probability should be estimated or \code{NULL}. If \code{NULL} the survival is estimated at all points in time at which an event occured.
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
I-Ming Chang, Rebecca Gelman, and Marcello Pagano. "Corrected Group Prognostic Curves and Summary Statistics". In: Journal of Chronic Diseases 35 (1982), pp. 669-674.

Robert W. Makuch. "Adjusted Survival Curve Estimation Using Covariates". In: Journal of Chronic Diseases 35.6 (1982), pp. 437-443.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
