\name{plot.adjusted_rmst}
\alias{plot.adjusted_rmst}

\title{
Plot Method for \code{adjusted_rmst} Objects
}
\description{
Plots the adjusted restricted mean survival times next to each other.
}
\usage{
\method{plot}{adjusted_rmst}(x, draw_ci=T, color=T, point_size=2,
     ci_size=1, ci_width=0.7, xlab="",
     ylab="Adjusted RMST", ...)
}

\arguments{
  \item{x}{
An object of class \code{adjusted_rmst} created by the \code{adjusted_rmst} function.
  }
  \item{draw_ci}{
Whether to draw the confidence intervals or not. Ignored if \code{use_boot=FALSE} in original call to the \code{adjusted_rmst} function. See \code{?adjusted_rmst} to get more information on confidence interval calculations.
  }
  \item{color}{
Whether to color the points (and if applicable confidence intervals) differently for each level of the grouping variable.
  }
  \item{point_size}{
The size of the points which will be drawn.
  }
  \item{ci_size}{
The size of the confidence interval bars. Ignored if \code{draw_ci=FALSE}.
  }
  \item{ci_width}{
The width of the confidence interval bars. Ignored if \code{draw_ci=FALSE}.
  }
  \item{xlab}{
The label of the X-Axis. Empty by default.
  }
  \item{ylab}{
The label of the Y-Axis. Defaults to "Adjusted RMST".
  }
  \item{title}{
The title of the plot. Set to \code{NULL} to use no title.
  }
  \item{...}{
Currently not used.
  }
}
\value{
Returns a \code{ggplot2} object.
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
