\name{plot.curve_test}
\alias{plot.curve_test}

\title{
Plot Method for \code{curve_test} Objects
}
\description{
Produces either a spaghetti-plot of the bootstrapped difference curves (\code{type="curves"}) or a density plot of the shifted bootstrap distribution of the difference curve integrals (\code{type="integral"}).
}
\usage{
\method{plot}{curve_test}(x, type="curves", ...)
}

\arguments{
  \item{x}{
An object of class \code{curve_test} created by the \code{test_curve_equality} function.
  }
  \item{type}{
Either \code{"curves"} or \code{"integral"}, specifying what should be plotted.
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
