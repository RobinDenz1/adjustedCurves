\name{print.adjusted_rmst}
\alias{print.adjusted_rmst}

\title{
Print Method for \code{adjusted_rmst} Objects
}
\description{
Prints some important parts of the output object.
}
\usage{
\method{print}{adjusted_rmst}(x, digits=5, ...)
}

\arguments{
  \item{x}{
An object of class \code{adjusted_rmst} created by the \code{adjusted_rmst} function.
  }
  \item{digits}{
How many digits to use when rounding the results.
  }
  \item{...}{
Currently not used.
  }
}
\value{
Silently returns the \code{data.frame} which can be seen when calling the function.
}
\author{
Robin Denz
}
\seealso{
\code{\link{adjusted_rmst}}, \code{\link{adjustedsurv}}
}
\examples{
# See ?adjusted_rmst
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
