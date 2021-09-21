\name{plot.adjustedsurv}
\alias{plot.adjustedsurv}

\title{
Plot Confounder-Adjusted Survival Curves
}
\description{
A function to graphically display confounder-adjusted survival curves which where previously estimated using the \code{adjustedsurv} function. The user can customize the plot using a variety of options. Internally it uses the \code{ggplot2} package, so additional not implemented features can be added using the standard \code{ggplot2} syntax. This function also includes the option to use isotonic regression on the survival curves, which is of benefit if the estimated curves are not monotone.
}
\usage{
\method{plot}{adjustedsurv}(x, draw_ci=F, max_t=Inf,
     iso_reg=F, force_bounds=F, use_boot=F,
     color=T, linetype=F, facet=F,
     line_size=1, xlab="Time",
     ylab="Adjusted Survival Probability",
     title=NULL, legend.title="Group",
     legend.position="right",
     ylim=NULL, custom_colors=NULL,
     custom_linetypes=NULL,
     ci_draw_alpha=0.4, steps=T, ...)
}

\arguments{
  \item{x}{
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
  \item{legend.title}{
A character string to be used as the title of the legend. Set to \code{NULL} if no legend should be included.
  }
  \item{legend.position}{
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
  \item{steps}{
A logical variable indicating whether the survival curves should be plotted as a step function or using straight lines. Straight lines should not be used with a simple Aalen-Joahnsen estimator. It is recommended to only use straight lines when a sufficiently fine grid of time points was used in the estimation step.
  }
  \item{...}{
Currently not used.
  }
}
\details{
When using certain methods there is no guarantee that the resulting estimated survival curves are monotonically decreasing. This is unfortunate since we know that it has to be the case. Isotonic regression can be used to fix this problem by ensuring that the survival curves are actually monotonically decreasing everywhere, while also being as close to the observations as possible. Westling et al. (2020) showed mathematically that this does not add any systematic bias to the estimats. More information on the method can be found in Robertson et al. (1988). This adjustment can be done using this function by setting \code{iso_reg} to \code{TRUE}.

Similarly, some methods can produce estimates that lie outside the theoretical 0 and 1 bounds of probability. By setting \code{force_bounds} to \code{TRUE} these estimates are manually set to either 0 or 1 (whichever is closer).
}
\value{
Returns a \code{ggplot2} object.
}
\references{
Ted Westling, Mark J. van der Laan, and Marco Carone. "Correcting an Estimator of a Multivariate Monotone Function with Isotonic Regression". In: Electronic Journal of Statistics 14 (2020), pp. 3032-3069.

Robertson, T., Wright, F. T. and Dykstra, R. L. "Order Restricted Statistical Inference". John Wiley & Sons, (1988)
}
\author{
Robin Denz
}

\seealso{
\code{\link{adjustedsurv}}, \code{\link[ggplot2]{ggplot}}, \code{\link[pammtools]{pammtools}}, \code{\link[stats]{isoreg}}
}
\examples{
\dontrun{
library(survival)
library(ggplot2)

# simulate some data as example
sim_dat <- sim_confounded_surv(n=500, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# estimate a cox-regression for the outcome
cox_mod <- coxph(Surv(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6 + group,
                 data=sim_dat, x=TRUE)


# use it to calculate adjusted survival curves with bootstrapping
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="direct",
                        outcome_model=cox_mod,
                        conf_int=TRUE,
                        bootstrap=TRUE,
                        n_boot=500)

# plot the curves with default values
plot(adjsurv)

# plot after applying isotonic regression
plot(adjsurv, iso_reg=TRUE)

# plot with confidence intervals estimated using asymptotic variances
plot(adjsurv, draw_ci=TRUE)

# plot with confidence intervals estimated using bootstrapping
plot(adjsurv, draw_ci=TRUE, use_boot=TRUE)

# plot with different linetypes only
plot(adjsurv, linetype=TRUE, color=FALSE, facet=FALSE)

# plot with different facets only
plot(adjsurv, linetype=FALSE, color=FALSE, facet=TRUE)

# plot with different linetypes and different colors
plot(adjsurv, linetype=TRUE, color=TRUE, facet=FALSE)

# plot with some custom characteristics
plot(adjsurv, legend.position="bottom", linetype=TRUE,
     custom_colors=c("green", "blue"), legend.title="Custom",
     title="Custom Plot", draw_ci=TRUE, linesize=0.5)

# plot with not directly included ggplot2 functionality
plot(adjsurv) + theme_classic()
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
