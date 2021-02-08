\name{adjustedsurv}
\alias{adjustedsurv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Confounder-Adjusted Survival Curves
}
\description{
This is the main function of this R-Package. It allows the user to calculate confounder-adjusted survival curves using any of the implemented methods. This function is basically a wrapper for each \code{surv_method_} function with some additional functionality, such as confidence intervals and bootstrapping. Depending on the specified method there might be additional required arguments. Documentation on that can be found in the respective \code{surv_method_} function. See details and examples.
}
\usage{
adjustedsurv(data, variable, ev_time, event, method, sd=T,
             times=NULL, alpha=0.05, bootstrap=F,
             n_boot=500, na.rm=F,...)
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
  \item{method}{
A character string specifying the adjustment method to use. See details.
  }
  \item{sd}{
A logical variable, indicating whether the standard deviation of the survival probabilities should be calculated. Needed to calculate confidence intervals.
  }
  \item{times}{
A numeric vector of time points at which the survival probability should be estimated or \code{NULL}. If \code{NULL} the survival is estimated at all points in time at which an event occured. Only relevant for some methods, ignored otherwise. See details.
  }
  \item{alpha}{
A number specifying the alpha error used in calculation of confidence intervals. Ignored if \code{sd=FALSE}.
  }
  \item{bootstrap}{
A logical variable indicating whether bootstrapping should be performed or not. See details.
  }
  \item{n_boot}{
Number of bootstrap replications to perform. Ignored if \code{bootstrap} is \code{FALSE}.
  }
  \item{na.rm}{
Whether to remove missing values or not. Only used in calculation of bootstrap statistics and methods using pseudo-values.
  }
  \item{...}{
Further arguments passed to the respective \code{surv_method}. See details.
  }
}
\details{

}
\value{
Returns an \code{adjustedsurv} object containing the following objects:

\item{adjsurv}{A \code{data.frame} of estimated adjusted survival probabilities for some points in time for each level of \code{variable}. Depending on the arguments used also includes standard deviations and confidence intervals at points in time.
}
\item{n}{The sample size
}
\item{method}{
The method used to adjust the survival curves.
}
\item{call}{
The original function call.
}

When the argument \code{bootstrap} is set to \code{TRUE} it additionally contains the following objects:

\item{boot_data}{The adjusted survival curves estimated in each bootstrap sample.
}
\item{boot_data_same_t}{The value of the adjusted survival probabilities at the same points in time as the original survival curve for each bootstrap sample.
}
\item{boot_adjsurv}{The mean adjusted survival curves of all bootstrap samples and corresponding standard deviations and percentile confidence intervals.
}

}

\references{
Hier meine eigenen Paper
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
