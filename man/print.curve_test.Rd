\name{print.curve_test}
\alias{print.curve_test}

\title{
Print Method for \code{curve_test} Objects
}
\description{
Prints some important parts of the output object.
}
\usage{
\method{print}{curve_test}(x, digits=4, ...)
}

\arguments{
  \item{x}{
An object of class \code{curve_test} created by the \code{adjusted_curve_diff} function.
  }
  \item{digits}{
How many digits to use when rounding the results.
  }
  \item{...}{
Currently not used.
  }
}
\value{
Silently returns the \code{data.frame} which can be seen when calling the function. \code{ABC} is an abbreviation for "area between the curves".
}
\author{
Robin Denz
}
\seealso{
\code{\link{adjusted_curve_diff}}, \code{\link{adjustedsurv}}, \code{\link{adjustedcif}}
}
\examples{
# See ?adjusted_curve_diff
}
