\name{plot.adjustedcif}
\alias{plot.adjustedcif}

\title{
Plot Confounder-Adjusted Cumulative Incidence Functions
}
\description{
A function to graphically display confounder-adjusted cumulative incidence functions which where previously estimated using the \code{\link[adjustedCurves]{adjustedcif}} function. The user can customize the plot using a variety of options. Internally it uses the \pkg{ggplot2} package, so additional not implemented features can be added using the standard \pkg{ggplot2} syntax. This function also includes the option to use isotonic regression on the CIFs, which is of benefit if the estimated curves are not monotone.
}
\usage{
\method{plot}{adjustedcif}(x, conf_int=FALSE, max_t=Inf,
     iso_reg=FALSE, force_bounds=FALSE,
     use_boot=FALSE, color=TRUE,
     linetype=FALSE, facet=FALSE,
     line_size=1, line_alpha=1, xlab="Time",
     ylab="Adjusted Cumulative Incidence",
     title=NULL, subtitle=NULL, legend.title="Group",
     legend.position="right",
     gg_theme=ggplot2::theme_classic(),
     ylim=NULL, custom_colors=NULL,
     custom_linetypes=NULL,
     single_color=NULL, single_linetype=NULL,
     conf_int_alpha=0.4, steps=TRUE,
     censoring_ind="none",
     censoring_ind_size=0.5,
     censoring_ind_alpha=1,
     censoring_ind_shape=17,
     censoring_ind_width=NULL,
     ...)
}

\arguments{
  \item{x}{
An \code{adjustedcif} object created using the \code{\link[adjustedCurves]{adjustedcif}} function.
  }
  \item{conf_int}{
A logical variable indicating whether the confidence intervals should be drawn.
  }
  \item{max_t}{
A number indicating the latest event time which is to be plotted.
  }
  \item{iso_reg}{
A logical variable indicating whether the estimates should be monotonized using isotonic regression. See details.
  }
  \item{force_bounds}{
A logical variable indicating whether the 0 and 1 bounds of the CIFs should be forced in the plot. See details.
  }
  \item{use_boot}{
A logical variable denoting whether the bootstrapped estimates should be used for the curves and their confidence intervals. Can only be used if they were calculated. See \code{\link[adjustedCurves]{adjustedcif}}.
  }
  \item{color}{
A logical variable indicating whether the curves should be colored differently. The \code{custom_colors} argument can be used to directly specify which colors to use. Alternatively the \code{single_color} argument can be used if everything should have the same color.
  }
  \item{linetype}{
A logical variable indicating whether the curves should have different linetypes. The \code{custom_linetypes} argument can be used to directly specify which linetypes to use. Alternatively the \code{single_linetype} argument can be used if all curves should have the same linetype.
  }
  \item{facet}{
A logical variable indicating whether the curves should be in different facets.
  }
  \item{line_size}{
A number controlling the thickness of the curves.
  }
  \item{line_alpha}{
A number controlling the transparency level of the curves.
  }
  \item{xlab}{
A character string to be used as the X-Axis label of the plot.
  }
  \item{ylab}{
A character string to be used as the Y-Axis label of the plot.
  }
  \item{title}{
A character string to be used as the title of the plot. Set to \code{NULL} (default) if no title should be used.
  }
  \item{subtitle}{
A character string to be used as the subtitle of the plot. Set to \code{NULL} (default) if no subtitle should be used.
  }
  \item{legend.title}{
A character string to be used as the title of the legend. Set to \code{NULL} if no legend should be included.
  }
  \item{legend.position}{
A character string specifying the position of the legend. Ignored if \code{legend_title=NULL}.
  }
  \item{gg_theme}{
A \code{ggplot2} theme object which will be used for the plot.
  }
  \item{ylim}{
A numeric vector of length two, specifying the limits of the Y-Axis. Set to \code{NULL} to use the \code{ggplot2} default values.
  }
  \item{custom_colors}{
A (named) vector to specify the colors of each CIF and possibly its confidence region. Set to \code{NULL} to use the \code{ggplot2} default values. Ignored if \code{color=FALSE}.
  }
  \item{custom_linetypes}{
A (named) vector to specify the linetype of each CIF. Set to \code{NULL} to use the \code{ggplot2} default values. Ignored if \code{linetype=FALSE}.
  }
  \item{single_color}{
A single color to use for every curve, irrespective of group status. If \code{color} is specified as well this argument will override it, but also generate a warning. Set to \code{NULL} (default) to ignore this argument.
  }
  \item{single_linetype}{
A single linetype to use for every curve, irrespective of group status. If \code{linetype} is specified as well this argument will override it, but also generate a warning. Set to \code{NULL} (default) to ignore this argument.
  }
  \item{conf_int_alpha}{
A number indicating the level of transparency that should be used when drawing the confidence regions.
  }
  \item{steps}{
A logical variable indicating whether the CIFs should be plotted as a step function or using straight lines. Straight lines should not be used with a simple Aalen-Joahnsen estimator. It is recommended to only use straight lines when a sufficiently fine grid of time points was used in the estimation step.
  }
  \item{censoring_ind}{
What kind of indicator to plot for censored observations on the CIFs. Must be one of \code{"none"} (plotting no indicators at all, the default), \code{"lines"} (plotting small vertical lines) and \code{"points"} (plotting points). Those will be affected by \code{linetype} and \code{color} as well. Observations who failed due to a competing event are not considered as censored here.
  }
  \item{censoring_ind_size}{
A numeric value specifying the size of the censoring indicators. Ignored if \code{censoring_ind="none"}.
  }
  \item{censoring_ind_alpha}{
A numeric value specifying the alpha level of the censoring indicators. Ignored if \code{censoring_ind="none"}.
  }
  \item{censoring_ind_shape}{
A numeric value specifying the shape of the censoring indicators when using \code{censoring_ind="points"}. Ignored otherwise. For available shapes see \code{?geom_point}.
  }
  \item{censoring_ind_width}{
A numeric value specifying the width of the censoring indicators. Ignored unless \code{censoring_ind="lines"}. By default (\code{censoring_ind_width=NULL}) the width of the censoring indicators is equal to 5 percent of the plot height.
  }
  \item{...}{
Currently not used.
  }
}
\details{
When using certain methods there is no guarantee that the resulting estimated CIFs are monotonically increasing. This is unfortunate since we know that it has to be the case. Isotonic regression can be used to fix this problem by ensuring that the CIFs are actually monotonically increasing everywhere, while also being as close to the observations as possible. Westling et al. (2020) showed mathematically that this usually does not add any systematic bias to the estimates. More information on the method can be found in Robertson et al. (1988). This adjustment can be done using this function by setting \code{iso_reg} to \code{TRUE}.

Similarly, some methods can produce estimates that lie outside the theoretical 0 and 1 bounds of probability. By setting \code{force_bounds} to \code{TRUE} these estimates are manually set to either 0 or 1 (whichever is closer).
}
\value{
Returns a \code{ggplot2} object.
}
\references{
Ted Westling, Mark J. van der Laan, and Marco Carone (2020). "Correcting an Estimator of a Multivariate Monotone Function with Isotonic Regression". In: Electronic Journal of Statistics 14, pp. 3032-3069

Tim Robertson, F. T. Wright, and R. L. Dykstra (1988). Order Restricted Statistical Inference. Hoboken: John Wiley & Sons
}
\author{
Robin Denz
}

\seealso{
\code{\link{adjustedcif}}, \code{\link[ggplot2]{ggplot}}, \code{\link[pammtools]{geom_stepribbon}}, \code{\link[stats]{isoreg}}
}
\examples{
library(adjustedCurves)
library(survival)

if (requireNamespace("riskRegression") & requireNamespace("prodlim") &
    requireNamespace("ggplot2")) {

library(riskRegression)
library(prodlim)
library(ggplot2)

set.seed(42)

# simulate some data as example
sim_dat <- sim_confounded_crisk(n=50)
sim_dat$group <- as.factor(sim_dat$group)

# calculate a Cause-Specific-Cox model
cox_mod <- CSC(Hist(time, event) ~ x1 + x3 + x5 + group,
               data=sim_dat)

# use it to calculate adjusted CIFs with bootstrapping (for cause = 1)
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="direct",
                      outcome_model=cox_mod,
                      conf_int=TRUE,
                      bootstrap=TRUE,
                      n_boot=15, # should be much bigger in reality
                      cause=1)

# plot the curves with default values
plot(adjcif)

# plot after applying isotonic regression
plot(adjcif, iso_reg=TRUE)

# plot with confidence intervals estimated using asymptotic variances
plot(adjcif, conf_int=TRUE)

# plot with confidence intervals estimated using bootstrapping
plot(adjcif, conf_int=TRUE, use_boot=TRUE)

# plot with different linetypes only
plot(adjcif, linetype=TRUE, color=FALSE, facet=FALSE)

# plot with different facets only
plot(adjcif, linetype=FALSE, color=FALSE, facet=TRUE)

# plot with different linetypes and different colors
plot(adjcif, linetype=TRUE, color=TRUE, facet=FALSE)

# plot with some custom characteristics
plot(adjcif, legend.position="bottom", linetype=TRUE,
     custom_colors=c("green", "blue"), legend.title="Custom",
     title="Custom Plot", conf_int=TRUE, linesize=0.5)

# adding further ggplot2 elements
plot(adjcif) + theme_bw()
}
}
