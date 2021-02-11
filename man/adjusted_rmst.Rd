\name{adjusted_rmst}
\alias{adjusted_rmst}

\title{
Calculate Confounder-Adjusted Restricted Mean Survival Times
}
\description{
This function can be utilized to calculate the confounder-adjusted restricted mean survival time, given previously estimated adjusted survival curves.
}
\usage{
adjusted_rmst(adjsurv, from=0, to=Inf, use_boot=F, conf_level=0.95)
}

\arguments{
  \item{adjsurv}{
An \code{adjustedsurv} object created using the \code{adjustedsurv} function.
  }
  \item{from}{
A number specifying the left side of the time interval of interest. See details.
  }
  \item{to}{
A number specifying the right side of the time interval of interest. See details.
  }
  \item{use_boot}{
Whether bootstrap estimates should be used to calculate the standard deviation of the RMST estimates.
  }
  \item{conf_level}{
A number specifying the confidence level of asymptotic and/or bootstrap confidence intervals.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a \code{adjusted_rmst} object containing the following object:

\item{areas}{A named vector containing the actual areas under the survival curves.
}
\item{rmsts}{A named vector containing the adjusted restricted mean survival time estimates.
}
\item{from}{\code{from} argument in function call.
}
\item{to}{\code{to} argument in function call.
}

If \code{use_boot} is set to \code{TRUE} and \code{bootstrap} was also set to \code{TRUE} in the original \code{adjustedsurv} function call, the following objects are also added to the output:

\item{booted_areas}{A \code{datal.frame} containing the areas under the survival curves for each bootstrap sample.
}
\item{booted_rmsts}{A \code{datal.frame} containing the adjusted restricted mean survival time estimates for each bootstrap sample.
}
\item{areas_sd}{A named vector containing the bootstrapped standard deviations of the areas under the survival curves.
}
\item{rmsts_sd}{A named vector containing the bootstrapped standard deviations of the restricted mean survival times.
}

}
\references{
Sarah C. Conner, Lisa M. Sullivan, Emelia J. Benjamin, Michael P. LaValley, Sandro Galea, and Ludovic Trinquart. "Adjusted Restricted Mean Survival Times in Observational Studies". In: Statistics in Medicine 38 (2019), pp. 3832-3860.

Patrick Royston and Mahesh K. B. Parmar. "Restricted Mean Survival Time: An Alternative to the Hazard Ratio for the Design and Analysis of Randomized Trials with a Time-To-Event Outcome". In: BMC Medical Research Methodology 13.152 (2013).
}
\author{
Robin Denz
}

\seealso{
\code{\link{print.adjusted_rmst}}, \code{\link{adjustedsurv}}
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
