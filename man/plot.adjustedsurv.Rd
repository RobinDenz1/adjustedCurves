\name{plot.adjustedsurv}
\alias{plot.adjustedsurv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Confounder-Adjusted Survival Curves
}
\description{
A function to graphically display confounder-adjusted survival curves which where previously estimated using the \code{adjustedsurv} function. The user can customize the plot using a variety of options. Internally it uses the \code{ggplot2} package, so additional not implemented features can be added using the standard \code{ggplot2} syntax. This function also includes the option to use isotonic regression on the survival curves, which is of benefit if the estimated curves are not monotone.
}
\usage{
plot.adjustedsurv(adjsurv, draw_ci=T, max_t=Inf,
                  iso_reg=F, force_bounds=F, use_boot=F,
                  color=T, linetype=F, facet=F,
                  line_size=1, xlab="Time",
                  ylab="Adjusted Survival Probability",
                  title=NULL, legend_title="Group",
                  legend_position="right",
                  ylim=NULL, custom_colors=NULL,
                  custom_linetypes=NULL,
                  ci_draw_alpha=0.4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adjsurv}{
An \code{adjustedsurv} object created using the \code{adjustedsurv} function.
  }
  \item{draw_ci}{
A logical variable indicating whether the confidence intervals should be drawn.
  }
  \item{max_t}{
A number indicating the latest survival time which is to be plotted.
  }
  \item{iso_reg}{
A logical variable indicating whether the estimates should be monotonized using isotonic regression. See details.
  }
  \item{force_bounds}{
A logical variable indicating whether the 0 and 1 bounds of the survival probabilities should be forced in the plot. See details.
  }
  \item{use_boot}{
A logical variable denoting whether the bootstrapped estimates should be used for the curves and their confidence intervals. Can only be used if they were calculated. See \code{adjustedsurv}.
  }
  \item{color}{
A logical variable indicating whether the curves should be colored differently.
  }
  \item{linetype}{
A logical variable indicating whether the curves should have different linetypes.
  }
  \item{facet}{
A logical variable indicating whether the curves should be in different facets.
  }
  \item{line_size}{
A number controlling the thickness of the survival curves.
  }
  \item{xlab}{
A character string to be used as the X-Axis label of the plot.
  }
  \item{ylab}{
A character string to be used as the Y-Axis label of the plot.
  }
  \item{title}{
A character string to be used as the title of the plot. Set to \code{NULL} if no title should be used.
  }
  \item{legend_title}{
A character string to be used as the title of the legend. Set to \code{NULL} if no legend should be included.
  }
  \item{legend_position}{
A character string specifying the position of the legend. Ignored if \code{legend_title=NULL}.
  }
  \item{ylim}{
A numeric vector of length two, specifying the limits of the Y-Axis. Set to \code{NULL} to use the \code{ggplot2} default values.
  }
  \item{custom_colors}{
A (named) vector to specify the colors of each adjusted survival curve and possibly its confidence region. Set to \code{NULL} to use the \code{ggplot2} default values.
  }
  \item{custom_linetypes}{
A (named) vector to specify the linetype of each adjusted survival curve. Set to \code{NULL} to use the \code{ggplot2} default values.
  }
  \item{ci_draw_alpha}{
A number indicating the level of transparency that should be used when drawing the confidence regions.
  }
}
\details{

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
\code{\link{adjustedsurv}}, \code{\link{ggplot2}}, \code{\link{pammtools}}
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
