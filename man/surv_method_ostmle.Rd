\name{surv_method_ostmle}
\alias{surv_method_ostmle}

\title{
One-Step Targeted Maximum Likelihood Estimation Survival Curves
}
\description{
This method utilizes the one-step targeted maximum likelihood methodology to estimate adjusted survival curves. It is very similar to \code{surv_method_tmle}, with some difference in the estimation techniques. Similar to the AIPTW methods it posseses the doubly-robust property and can also be used to adjust for dependent censoring. It is however important to note that this method should only be used when the time variable is truly discrete and only a few distinct points in time are observed. See details for more information.
}
\usage{
surv_method_ostmle(data, variable, ev_time, event, conf_int,
                   conf_level=0.95, times, t_max=NULL,
                   adjust_vars=NULL, SL.ftime=NULL, SL.ctime=NULL,
                   SL.trt=NULL, epsilon=1, max_num_iteration=100,
                   psi_moss_method="l2", tmle_tolerance=NULL,
                   gtol=1e-3)
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
  \item{t_max}{
The furthest point in time to consider in the survival curve estimation. See \code{MOSS}.
}
  \item{adjust_vars}{
A character vector of column names specifying variables to be used when modeling the outcome, treatment and censoring mechanism. See details and examples.
}
  \item{SL.ftime}{
A character vector or list specification to be passed to the \code{SL.library} option in the call to \code{SuperLearner} for the outcome regression. See\code{?SuperLearner} for more information on how to specify valid \code{SuperLearner} libraries. It is expected that the wrappers used in the library will play nicely with the input variables, which will be called "trt", \code{names(adjust_vars)}, and "t".
  }
  \item{SL.ctime}{
A character vector or list specification to be passed to the \code{SL.library} argument in the call to \code{SuperLearner} for the estimate of the conditional hazard for censoring. It is expected that the wrappers used in the library will play nicely with the input variables, which will be called "trt" and \code{names(adjust_vars)}.
  }
  \item{SL.trt}{
A character vector or list specification to be passed to the \code{SL.library} argument in the call to \code{SuperLearner} for the estimate of the conditional probability oftreatment. It is expected that the wrappers used in the library will play nicely with the input variables, which will be \code{names(adjust_vars)}.
  }
  \item{epsilon}{

}
  \item{max_num_iteration}{

}
  \item{psi_moss_method}{

}
  \item{tmle_tolerance}{

}
  \item{gtol}{

}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a \code{data.frame} object containing the estimated adjusted survival probabilities for some points in time for each level of \code{variable}.
}
\references{
Weixin Cai and Mark J. van der Laan. "One-Step Targeted Maximum Likelihood Estimation for Time-To-Event Outcomes". In: Biometrics 76 (2020), pp. 722-733.
}
\author{
The wrapper function was written by Robin Denz, the \code{MOSS} package (which this wrapper is based on) was written by Wilson Cai. See \code{?MOSS} for more details.
}

\seealso{
\code{\link[survtmle]{survtmle}}, \code{\link[MOSS]{MOSS}}, \code{\link[SuperLearner]{SuperLearner}}
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
