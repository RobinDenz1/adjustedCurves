\name{surv_method_aiptw_pseudo}
\alias{surv_method_aiptw_pseudo}

\title{
Augmented Inverse Probability of Treatment Weighted Survival Curves using Pseudo-Values
}
\description{
Function to calculate augmented inverse probability of treatment weighted survival curves using pseudo-values. Instead of only modeling the outcome mechanism or the treatment assignment mechanism, both kind of models are required to use this method. If either of those models are correctly specified, unbiased estimates will be obtained. Both the outcome and treatment assignment part are modeled using pseudo-values. See details for more information.
}
\usage{
surv_method_aiptw_pseudo(data, variable, ev_time, event, conf_int,
                         conf_level=0.95, times, outcome_vars,
                         treatment_model, type_time="factor",
                         spline_df=10, na.rm=F)
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
  \item{outcome_vars}{
A character vector of column names specifying variables to be used when modeling the outcome mechanism using \code{geese}. See details and examples.
}
  \item{treatment_model}{
Must be a \code{glm} model object with \code{variable} as response variable. See details and examples.
}
  \item{type_time}{
A character string specifying how the time should be modeled. See details.
}
  \item{spline_df}{
The number of degrees of freedom used for the natural-spline or B-spline function. Ignored if \code{type_time="factor"}.
}
  \item{na.rm}{
A logical variable, indicating whether rows with missing Pseudo-Values should be removed.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a \code{data.frame} object containing the estimated adjusted survival probabilities for some points in time for each level of \code{variable}.
}
\references{

Jixian Wang. "A Simple, Doubly Robust, Efficient Estimator for Survival Functions Using Pseudo Observations". In: Pharmaceutical Statistics 17.38-48 (2018).

James M. Robins and Andrea Rotnitzky. "Recovery of Information and Adjustment for Dependent Censoring Using Surrogate Markers". In: AIDS Epidemiology: Methodological Issues. Ed. by Nicholas P. Jewell, Klaus Dietz, and Vernon T. Farewell. New York: Springer Science + Business Media, 1992, pp. 297-331.

Per K. Andersen, Elisavet Syriopoulou, and Erik T. Parner. "Causal Inference in Survival Analysis using Pseudo-Observations". In: Statistics in Medicine 36 (2017), pp. 2669-2681.

Per K. Andersen and Maja Pohar Perme. "Pseudo-Observations in Survival Analysis". In: Statistical Methods in Medical Research 19 (2010), pp. 71-99.

Aris Perperoglou, Willi Sauerbrei, Michal Abrahamowicz, and Matthias Schmid. "A Review of Spline Function Procedures in R". In: BMC Medical Research Methodology 19 (2019), pp. 1-16
}
\author{
Robin Denz
}

\seealso{
\code{\link[geepack]{geese}}, \code{\link[prodlim]{jackknife}}, \code{\link[splines]{ns}}, \code{\link[splines]{bs}}
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
