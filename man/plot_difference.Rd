\name{plot_difference}
\alias{plot_difference}

\title{
Plot the Difference of Two Adjusted Survival Curves or CIFs
}
\description{
A function to graphically display the difference between two confounder-adjusted survival curves which where previously estimated using the \code{\link[adjustedCurves]{adjustedsurv}} function or between two confounder-adjusted CIFs which where previously estimated using the \code{\link[adjustedCurves]{adjustedcif}} function. The user can customize the plot using a variety of options. Internally it uses the \pkg{ggplot2} package, so additional not implemented features can be added using the standard \code{ggplot2} syntax.
}
\usage{
plot_difference(x, group_1=NULL, group_2=NULL,
                conf_int=FALSE, conf_level=0.95,
                type="steps", max_t=Inf,
                size=0.7, color="black", linetype="solid",
                alpha=1, conf_int_alpha=0.4,
                points_ci_size=NULL, points_ci_width=NULL,
                xlab="Time", ylab=NULL, title=NULL,
                subtitle=NULL, gg_theme=ggplot2::theme_classic(),
                line_at_0=TRUE, line_at_0_size=0.7,
                line_at_0_color="grey", line_at_0_linetype="dashed",
                line_at_0_alpha=1,
                loess_smoother=FALSE, loess_span=0.75,
                loess_color=color, loess_size=size,
                loess_linetype="dashed", loess_alpha=alpha,
                ...)
}

\arguments{
  \item{x}{
An \code{adjustedsurv} object created using the \code{\link[adjustedCurves]{adjustedsurv}} function or an \code{adjustedcif} object created using the \code{\link[adjustedCurves]{adjustedcif}} function.
  }
  \item{group_1}{
A single character string specifying one of the levels of the \code{variable} used in the original \code{adjustedsurv} or \code{adjustedcif} function call. This group will be subtracted from. For example if \code{group_1="A"} and \code{group_2="B"} the plotted curve will correspond to the survival probability of \code{A} minus the survival probability of \code{B} over time.
  }
  \item{group_2}{
Also a single character string specifying one of the levels of \code{variable}. This corresponds to the right side of the difference equation. See argument \code{group_2}.
  }
  \item{conf_int}{
A logical variable indicating whether the confidence intervals should be drawn. This only works when \code{bootstrap=TRUE} was used in the original \code{adjustedsurv} or \code{adjustedcif} function call. It also currently does not work if multiply imputed data was used.
  }
  \item{conf_level}{
The confidence level that should be used when calculating the confidence intervals. Ignored if \code{conf_int=FALSE}.
  }
  \item{type}{
Must be one of \code{"steps"} (drawing the difference as a step function), \code{"lines"} (drawing the difference using linear interpolation), \code{"points"} (drawing points only) or \code{"none"} (drawing nothing, useful when only the smoothed difference is of interest). It defaults to \code{"steps"}.
  }
  \item{max_t}{
A number indicating the latest time to which the curve should be extended to.
  }
  \item{size}{
A number controlling the thickness of the difference curve.
  }
  \item{color}{
A string specifying the color of the difference curve.
  }
  \item{linetype}{
A string specifying the linetype of the difference curve.
  }
  \item{alpha}{
A number controlling the transparency level of the difference curves.
  }
  \item{conf_int_alpha}{
A number indicating the level of transparency that should be used when drawing the confidence regions.
  }
  \item{points_ci_size}{
Only used when \code{type="points"}. Controls the size of the error bars.
  }
  \item{points_ci_width}{
Only used when \code{type="points"}. Controls the width of the error bars.
  }
  \item{xlab}{
A character string to be used as the X-Axis label of the plot. Defaults to "Time".
  }
  \item{ylab}{
A character string to be used as the Y-Axis label of the plot. By default (\code{NULL}) uses the equation used to calculate the differences, based on the names supplied in \code{group_1} and \code{group_2}.
  }
  \item{title}{
A character string to be used as the title of the plot. Set to \code{NULL} if no title should be used.
  }
  \item{subtitle}{
A character string to be used as the subtitle of the plot. Set to \code{NULL} if no subtitle should be used.
  }
  \item{gg_theme}{
A \code{ggplot2} theme object which will be used for the plot.
  }
  \item{line_at_0}{
Whether to draw a horizontal line at y = 0 or not.
  }
  \item{line_at_0_size}{
The size of the line drawn at y = 0. Ignored if \code{line_at_0=FALSE}.
  }
  \item{line_at_0_color}{
The color of the line drawn at y = 0. Ignored if \code{line_at_0=FALSE}.
  }
  \item{line_at_0_linetype}{
The linetype of the line drawn at y = 0. Ignored if \code{line_at_0=FALSE}.
  }
  \item{line_at_0_alpha}{
The transparency level of the line drawn at y = 0. Ignored if \code{line_at_0=FALSE}.
  }
  \item{loess_smoother}{
Whether to draw a LOESS smoother through the difference curves.
  }
  \item{loess_span}{
The span of the LOESS smoother. Ignored if \code{loess_smoother=FALSE}. See \code{stat_smooth} in the \code{ggplot2} package, \code{method="loess"} for more details.
  }
  \item{loess_color}{
The color of the LOESS smoother line. Ignored if \code{loess_smoother=FALSE}.
  }
  \item{loess_size}{
The size of the LOESS smoother line. Ignored if \code{loess_smoother=FALSE}.
  }
  \item{loess_linetype}{
The linetype of the LOESS smoother line. Ignored if \code{loess_smoother=FALSE}.
  }
  \item{loess_alpha}{
The transparency level of the LOESS smoother line. Ignored if \code{loess_smoother=FALSE}.
  }
  \item{...}{
Currently not used.
  }
}
\details{
This function allows the easy creation of difference curves. The syntax is exactly the same for both adjusted survival curves and adjusted CIFs. It currently does not support plotting multiple difference curves at once, which could be useful when there are more than two treatment groups in \code{variable}. If the user is interested in this, we recommend calling this function multiple times with the desired comparisons and concatenating the individual plots into one plot afterwards using a suitable function such as \code{par} or \code{ggarrange}.
}
\value{
Returns a \code{ggplot2} object.
}
\author{
Robin Denz
}

\seealso{
\code{\link{adjustedsurv}}, \code{\link{adjustedcif}}, \code{\link[ggplot2]{ggplot}}, \code{\link[pammtools]{pammtools}}
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

# plot the difference with default values
plot_difference(adjsurv)

# plot with reversed differences
plot_difference(adjsurv, group_1="1", group_2="0")

# plot with confidence intervals
plot_difference(adjsurv, conf_int=TRUE)

# plot using lines instead
plot_difference(adjsurv, conf_int=TRUE, type="lines")

# plot using points instead
plot_difference(adjsurv, conf_int=TRUE, type="points")

# plot using an additional loess smoother
plot_difference(adjsurv, loess_smoother=TRUE)

# plot without the line at 0
plot_difference(adjsurv, line_at_0=FALSE)

# plot with some custom parameters
plot_difference(adjsurv, conf_int=TRUE, color="blue", linetype="dotted",
                alpha=0.8, line_at_0_size=1.1, line_at_0_color="red",
                loess_smoother=TRUE, loess_span=0.55)

}
}
