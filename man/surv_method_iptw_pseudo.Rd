\name{surv_method_iptw_pseudo}
\alias{surv_method_iptw_pseudo}

\title{
Inverse Probability of Treatment Weighted Survival Estimates using Pseudo-Values
}
\description{
This method works by modelling the treatment assignment mechanism. The individual observations are weighted using inverse probability of treatment weights or other balancing weights (see details) and pseudo-values. If the model is correctly specified, unbiased results will be obtained. The main difference to the \code{surv_method_iptw_cox} and \code{surv_method_iptw_km} functions is the use of pseudo-values. See details for more information.
}
\usage{
surv_method_iptw_pseudo(data, variable, ev_time, event, conf_int,
                        conf_level=0.95, times, treatment_model,
                        weight_method="ps", se_method="cochrane", na.rm=F, ...)
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
A logical variable, indicating whether the asymptotic variances and confidence intervals of the survival probabilities should be calculated.
  }
  \item{conf_level}{
A number specifying the confidence level of asymptotic and/or bootstrap confidence intervals.
  }
  \item{times}{
A numeric vector of time points at which the survival probability should be estimated or \code{NULL}. If \code{NULL} the survival is estimated at all points in time at which an event occured.
}
  \item{treatment_model}{
Must be either a model object with \code{variable} as response variable, a vector of weights or a formula which can be passed to \code{WeightIt}.
}
  \item{weight_method}{
Method used in \code{WeightIt} function call. Ignored if \code{treatment_model} is not a formula object.
}
  \item{se_method}{
One of \code{"miller"}, \code{"galloway"} and \code{"cochrane"}. Specifies which kind of standard error to calculate. Defaults to \code{"cochrane"}. See details.
}
  \item{na.rm}{
A logical variable, indicating whether rows with missing Pseudo-Values should be removed.
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
Per K. Andersen, Elisavet Syriopoulou, and Erik T. Parner. "Causal Inference in Survival Analysis using Pseudo-Observations". In: Statistics in Medicine 36 (2017), pp. 2669-2681.

Per K. Andersen and Maja Pohar Perme. "Pseudo-Observations in Survival Analysis". In: Statistical Methods in Medical Research 19 (2010), pp. 71-99.

Donald F. Gatz and Luther Smith. "The Standard Error of a Weighted Mean Concentration - I: Bootstrapping Vs Other Methods". In: Atmospheric Environment 29.11 (1995), pp. 1185-1193.

William G. Cochran. Sampling Techniques. Vol. 3. New York: Wiley, 1977.

J. N. Galloway, G. E. Likens, and M. E. Hawley. "Acid Precipitation: Natural Versus Anthropogenic Components". In: Science 226 (1984), pp. 829-831.

J. M. Miller. A Statistical Evaluation of the U.S. Precipitation Chemistry Network. Precipitation Scavenging (edited by Semonin R. G. and Beadle R. W.) pp. 639-659. Available as CONF 74100 from National Technical Information Service, U.S. Dept. of Commerce, Springfiel, VA. 1977.
}
\author{
Robin Denz
}

\seealso{
\code{\link{WeightIt}}, \code{\link{prodlim}}
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
