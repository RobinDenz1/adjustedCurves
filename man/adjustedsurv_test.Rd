\name{adjustedsurv_test}
\alias{adjustedsurv_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test for the Equality of Two Confounder-Adjusted Survival Curves
}
\description{
This function implements a modified version of the Pepe-Flemming (ZITAT) test for the equality of two adjusted survival curves. In particular, the Null-Hypothesis is that the integral of the difference of the two curves in a specified time interval is not different from zero. The input of this function should be a \code{adjustedsurv} object created using the \code{adjustedsurv} function, with \code{bootstap=TRUE} in the original function call.
}
\usage{
adjustedsurv_test(adjsurv, from=0, to=Inf)
}
%- maybe also 'usage' for other objects documented here.
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
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a \code{adjustedsurv_test} object containing the following object:

\item{diff_curves}{A \code{data.frame} containing the difference curves used for calculating the integrals.
}
\item{diff_intergals}{A numeric vector containing the integrals of the difference for the estimated adjusted survival curves.
}
\item{observed_diff_curve}{The curve of the difference between the two non-bootstrapped adjusted survival curves.
}
\item{obersved_diff_integral}{The integral of the curve in \code{observed_diff_curve}.
}
\item{p_value}{The P-Value for the modified Pepe-Fleming Test. See details.
}

}
\references{

}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{print.adjustedsurv_test}}, \code{\link{adjustedsurv}}
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
