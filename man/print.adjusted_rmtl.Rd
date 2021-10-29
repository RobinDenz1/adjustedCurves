\name{print.adjusted_rmtl}
\alias{print.adjusted_rmtl}

\title{
Print Method for \code{adjusted_rmtl} Objects
}
\description{
Prints some important parts of the output object.
}
\usage{
\method{print}{adjusted_rmtl}(x, digits=5, ...)
}

\arguments{
  \item{x}{
An object of class \code{adjusted_rmtl} created by the \code{adjusted_rmtl} function.
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
\code{\link{adjusted_rmtl}}, \code{\link{adjustedsurv}}, \code{\link{adjustedcif}}
}
\examples{
# See ?adjusted_rmtl
}
