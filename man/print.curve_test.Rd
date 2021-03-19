\name{print.curve_test}
\alias{print.curve_test}

\title{
Print Method for \code{curve_test} Objects
}
\description{
Prints some important parts of the output object.
}
\usage{
\method{print}{curve_test}(x, ...)
}

\arguments{
  \item{x}{
An object of class \code{curve_test} created by the \code{test_curve_equality} function.
  }
  \item{...}{
Currently not used.
  }
}

\value{
\code{NULL}
}

\author{
Robin Denz
}

\seealso{
\code{\link{test_curve_equality}}, \code{\link{adjustedsurv}}, \code{\link{adjusted_rmst}}
}
\examples{
# See ?test_curve_equality
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
