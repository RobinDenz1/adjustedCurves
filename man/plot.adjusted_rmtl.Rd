\name{plot.adjusted_rmtl}
\alias{plot.adjusted_rmtl}

\title{
Plot Method for \code{adjusted_rmtl} Objects
}
\description{
Plots the adjusted restricted mean time lost estimates next to each other.
}
\usage{
\method{plot}{adjusted_rmtl}(x, draw_ci=TRUE, color=TRUE, point_size=2,
     ci_size=1, ci_width=0.7, xlab="",
     ylab="Adjusted RMTL", title=NULL,
     gg_theme=ggplot2::theme_classic(), ...)
}

\arguments{
  \item{x}{
An object of class \code{adjusted_rmtl} created by the \code{adjusted_rmtl} function.
  }
  \item{draw_ci}{
Whether to draw the confidence intervals or not. Ignored if \code{use_boot=FALSE} in original call to the \code{adjusted_rmtl} function. See \code{?adjusted_rmtl} to get more information on confidence interval calculations.
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
The label of the Y-Axis. Defaults to "Adjusted RMTL".
  }
  \item{title}{
The title of the plot. Set to \code{NULL} to use no title.
  }
  \item{gg_theme}{
A \code{ggplot2} theme object which will be used for the plot.
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
\code{\link{adjusted_rmtl}}, \code{\link{adjustedsurv}}, \code{\link{adjustedcif}}
}
\examples{
# See ?adjusted_rmtl
}
