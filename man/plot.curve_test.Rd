\name{plot.curve_test}
\alias{plot.curve_test}

\title{
Plot Method for \code{curve_test} Objects
}
\description{
Produces either a spaghetti-plot of the bootstrapped difference curves (\code{type="curves"}) or a kernel-density plot of the shifted bootstrap distribution of the difference curve integrals (\code{type="integral"}).
}
\usage{
\method{plot}{curve_test}(x, type="curves", xlab=NULL, ylab=NULL, title=NULL, ...)
}

\arguments{
  \item{x}{
An object of class \code{curve_test} created by the \code{test_curve_equality} function.
  }
  \item{type}{
Either \code{"curves"} or \code{"integral"}, specifying what should be plotted.
  }
  \item{xlab}{
The label of the X-Axis. Set to \code{NULL} to use default label.
  }
  \item{ylab}{
The label of the Y-Axis. Set to \code{NULL} to use default label.
  }
  \item{title}{
The title of the plot. Set to \code{NULL} to use no title.
  }
  \item{...}{
Currently not used.
  }
}

\details{

When using \code{type="curves"} the black curve shows the observed curve of the difference. When using \code{type="integral"} the red line shows the observed integral of the curve of the difference.

Both graphics can be used to check if the assumptions of the test hold. The bootstrap-shifted distribution of the integral of the difference should approximately be normally distributed. If the kernel-density estimate shown with \code{type="integral"} is clearly not normally distributed, the estimated p-value might be wrong. Similarly, if the curves of the differences do not vary randomly around the black line when using \code{type="curves"}, the estimated p-value might be wrong. You could also try to rerun the \code{adjustedsurv} or \code{adjustedcif} function with a bigger number in \code{n_boot}.

}

\value{
Returns a \code{ggplot2} object.
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
