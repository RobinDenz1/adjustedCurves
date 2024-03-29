\name{plot_curve_diff}
\alias{plot_curve_diff}
\alias{plot_curve_ratio}

\title{
Plot the Difference Between or the Ratio of Two Adjusted Survival Curves or CIFs
}
\description{
A function to graphically display the difference or ratio between two confounder-adjusted survival curves which where previously estimated using the \code{\link[adjustedCurves]{adjustedsurv}} function or between two confounder-adjusted CIFs which where previously estimated using the \code{\link[adjustedCurves]{adjustedcif}} function. The user can customize the plot using a variety of options. Internally it uses the \pkg{ggplot2} package, so additional not implemented features can be added using the standard \code{ggplot2} syntax.
}
\usage{
plot_curve_diff(x, group_1=NULL, group_2=NULL,
                conf_int=FALSE, conf_level=0.95, type="steps",
                times=NULL, max_t=Inf, use_boot=FALSE,
                size=0.7, color="black", linetype="solid",
                alpha=1, conf_int_alpha=0.4,
                points_ci_size=NULL, points_ci_width=NULL,
                xlab="Time", ylab=NULL, title=NULL,
                subtitle=NULL, gg_theme=ggplot2::theme_classic(),
                line_at_ref=TRUE, line_at_ref_size=0.7,
                line_at_ref_color="grey", line_at_ref_linetype="dashed",
                line_at_ref_alpha=1,
                loess_smoother=FALSE, loess_span=0.75,
                loess_color=color, loess_size=size,
                loess_linetype="dashed", loess_alpha=alpha,
                test=NULL, integral_from=0, integral_to=NULL,
                p_value=FALSE, integral=FALSE,
                interval=FALSE, text_pos_x="left",
                text_pos_y="bottom", text_size=3.5,
                text_family="serif", text_fontface="italic",
                text_color="black", text_alpha=1,
                text_digits=3, text_format_p=TRUE,
                fill_area=FALSE, area_color="blue", area_alpha=0.4,
                fill_only_interval=TRUE,
                ...)

plot_curve_ratio(x, group_1=NULL, group_2=NULL, conf_int=FALSE,
                 conf_level=0.95, type="steps", times=NULL,
                 max_t=Inf, use_boot=FALSE, size=0.7, color="black",
                 linetype="solid", alpha=1,
                 conf_int_alpha=0.4, xlab="Time", ylab=NULL,
                 title=NULL, subtitle=NULL,
                 gg_theme=ggplot2::theme_classic(),
                 line_at_ref=TRUE, line_at_ref_size=0.7,
                 line_at_ref_color="grey",
                 line_at_ref_linetype="dashed",
                 line_at_ref_alpha=1, ...)
}

\arguments{
  \item{x}{
An \code{adjustedsurv} object created using the \code{\link[adjustedCurves]{adjustedsurv}} function or an \code{adjustedcif} object created using the \code{\link[adjustedCurves]{adjustedcif}} function.
  }
  \item{group_1}{
A single character string specifying one of the levels of the \code{variable} used in the original \code{adjustedsurv} or \code{adjustedcif} function call. This group will be subtracted from. For example if \code{group_1="A"} and \code{group_2="B"} the plotted curve will correspond to the survival probability (or CIF) of \code{A} minus the survival probability (or CIF) of \code{B} over time. If \code{NULL}, this will default to the first level of \code{variable}. Similarly if one plots ratios instead, the ratio would be calculated as \code{A / B}.
  }
  \item{group_2}{
Also a single character string specifying one of the levels of \code{variable}. This corresponds to the right side of the difference equation. See argument \code{group_2}. If \code{NULL}, this will default to the second level of \code{variable}.
  }
  \item{conf_int}{
A logical variable indicating whether the confidence intervals should be drawn. This only works when \code{conf_int=TRUE} or \code{bootstrap=TRUE} was used in the original \code{adjustedsurv} or \code{adjustedcif} function call.
  }
  \item{conf_level}{
The confidence level that should be used when calculating the confidence intervals. Ignored if \code{conf_int=FALSE}.
  }
  \item{type}{
Must be one of \code{"steps"} (drawing the difference/ratio as a step function), \code{"lines"} (drawing the difference/ratio using linear interpolation), \code{"points"} (drawing points only) or \code{"none"} (drawing nothing, useful when only the smoothed difference is of interest). It defaults to \code{"steps"}. For ratios, only the \code{"steps"} and \code{"lines"} options are available.
  }
  \item{times}{
An optional numeric vector of points in time at which the difference or ratio should be estimated. If \code{NULL} (default) the differences / ratios are estimated for the whole curve. This only affects the plot and has no effect on the \code{integral} or \code{p_value} if those are also specified.
  }
  \item{max_t}{
A number indicating the latest time to which the curve should be extended to.
  }
  \item{use_boot}{
Whether to use the bootstrapped estimates to calculate the confidence intervals or not. Can only be used if \code{bootstrap=TRUE} was used in the original \code{adjustedsurv} or \code{adjustedcif} function call. Ignored if \code{conf_int=FALSE}.
  }
  \item{size}{
A number controlling the thickness of the curve.
  }
  \item{color}{
A string specifying the color of the curve.
  }
  \item{linetype}{
A string specifying the linetype of the curve.
  }
  \item{alpha}{
A number controlling the transparency level of the curves.
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
A character string to be used as the Y-Axis label of the plot. By default (\code{NULL}) uses the equation used to calculate the differences / ratios, based on the names supplied in \code{group_1} and \code{group_2}.
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
  \item{line_at_ref}{
Whether to draw a horizontal line at y = 0 for differences or at y = 1 for ratios or not.
  }
  \item{line_at_ref_size}{
The size of the line drawn at the reference value. Ignored if \code{line_at_ref=FALSE}.
  }
  \item{line_at_ref_color}{
The color of the line drawn at the reference value. Ignored if \code{line_at_ref=FALSE}.
  }
  \item{line_at_ref_linetype}{
The linetype of the line drawn at the reference value. Ignored if \code{line_at_ref=FALSE}.
  }
  \item{line_at_ref_alpha}{
The transparency level of the line drawn at the reference value. Ignored if \code{line_at_ref=FALSE}.
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
  \item{test}{
An optional \code{curve_test} object created using the \code{adjusted_curve_test} function. If supplied it can be used to add a p-value and the integral statistic to the plot. Alternatively, the needed arguments below can be specified to obtain the values needed for the test. See below. Set to \code{NULL} (default) to ignore this.
  }
  \item{integral_from}{
A number specifying the left limit of the integral. When \code{p_value=TRUE} and \code{test=NULL}, this argument will be passed to the \code{from} argument in the \code{adjusted_curve_test} function to perform the test.
  }
  \item{integral_to}{
A number specifying the right limit of the integral. When \code{p_value=TRUE} and \code{test=NULL}, this argument will be passed to the \code{to} argument in the \code{adjusted_curve_test} function to perform the test.
  }
  \item{p_value}{
Whether to add a p-value to the plot or not. This requires either that the user supplies a previously created \code{curve_test} object to the \code{test} argument, or that the required arguments to call this function are supplied (at least \code{integral_to}). Either way it only works if \code{bootstrap=TRUE} was used in the original \code{adjustedsurv} or \code{adjustedcif} function call.
  }
  \item{integral}{
Whether to add the integral of the difference in the interval [\code{from}, \code{to}] to the plot or not. This requires either that the user supplies a previously created \code{curve_test} object to the \code{test} argument, or that the required arguments to call this function are supplied (at least \code{integral_to}).
  }
  \item{interval}{
Whether to add the interval in which the integral was calculated to the plot as well.
  }
  \item{text_pos_x}{
X position of the text. Can be either \code{"left"} (default), \code{"middle"}, \code{"right"} or a number specifying the exact position.
  }
  \item{text_pos_y}{
Y position of the text. Can be either \code{"bottom"} (default), \code{"middle"}, \code{"top"} or a number specifying the exact position.
  }
  \item{text_digits}{
The number of digits to which the p-value and the integral of the difference should be rounded to.
  }
  \item{text_size}{
The size of the text.
  }
  \item{text_family}{
The family of the text. Defaults to \code{"serif"}.
  }
  \item{text_fontface}{
The fontface of the text. Defaults to \code{"italic"}.
  }
  \item{text_color}{
The color of the text. Defaults to \code{"black"}.
  }
  \item{text_alpha}{
The transparency level of the text.
  }
  \item{text_format_p}{
Whether to format p-values smaller than 0.01 to < 0.01.
  }
  \item{fill_area}{
Whether to add color to the area between 0 and the difference.
  }
  \item{area_color}{
The color used to fill in the area between 0 and the difference when using \code{fill_area=TRUE}. Ignored otherwise.
  }
  \item{area_alpha}{
The transparency level used to fill in the area between 0 and the difference when using \code{fill_area=TRUE}. Ignored otherwise.
  }
  \item{fill_only_interval}{
Whether only the area corresponding to the interval defined by \code{integral_from} and \code{integral_to} should be filled. Only used when \code{fill_area=TRUE}.
  }
  \item{...}{
Currently not used.
  }
}
\details{
This function allows the easy creation of difference / ratio curves. The syntax is exactly the same for both adjusted survival curves and adjusted CIFs. Similarly, the syntax is the same for ratios and for difference curves, although not all options of the difference curve function are available for the ratio curve function. By default it calculates the estimates up to the last point where estimates for both the \code{group_1} curve and the \code{group_2} curve are available.

It currently does not support plotting multiple curves at once, which could be useful when there are more than two treatment groups in \code{variable}. If the user is interested in this, we recommend calling this function multiple times with the desired comparisons and concatenating the individual plots into one plot afterwards using a suitable function such as \code{par} or \code{ggarrange}.

More information on how the differences or ratios and their confidence intervals are calculated can be found in the documentation of the \code{\link{adjusted_curve_diff}} function. More information on how the overall p-value and the integral are calculated for differences can be found in the \code{\link{adjusted_curve_test}} function.
}
\value{
Returns a \code{ggplot2} object.
}
\references{
Michael Coory, Karen E. Lamb, and Michael Sorich (2014). "Risk-Difference Curves can be used to Communicate Time-Dependent Effects of Adjuvant Therapies for Early Stage Cancer". In: Journal of Clinical Epidemiology 67, pp. 966-972

Lihui Zhao, Lu Tian, Hajime Uno, Scott D. Solomon, Marc A. Pfeffer, Jerald S. Schindler, and L. J. Wei (2012). "Utilizing the Integrated Difference of Two Survival Functions to Quantify the Treatment Contrast for Designing, Monitoring and Analyzing a Comparative Clinical Study". In: Clinical Trials 9.5, pp. 570-577
}
\author{
Robin Denz
}

\seealso{
\code{\link{adjusted_curve_diff}}, \code{\link{adjusted_curve_test}}, \code{\link{adjustedsurv}}, \code{\link{adjustedcif}}, \code{\link[ggplot2]{ggplot}}, \code{\link[pammtools]{geom_stepribbon}}
}
\examples{
library(adjustedCurves)
library(survival)

if (requireNamespace("ggplot2") & requireNamespace("riskRegression")) {

library(ggplot2)

set.seed(42)

# simulate some data as example
sim_dat <- sim_confounded_surv(n=50, max_t=1.2)
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
                        n_boot=15) # should be much bigger in reality

# plot the difference with default values
plot_curve_diff(adjsurv)

# plot the ratio with default values
plot_curve_ratio(adjsurv)

# plot with reversed differences
plot_curve_diff(adjsurv, group_1="1", group_2="0")

# plot with confidence intervals
plot_curve_diff(adjsurv, conf_int=TRUE)
plot_curve_ratio(adjsurv, conf_int=TRUE)

# plot using lines instead
plot_curve_diff(adjsurv, conf_int=TRUE, type="lines")

# plot using points instead
plot_curve_diff(adjsurv, conf_int=TRUE, type="points")

# plot using an additional loess smoother
plot_curve_diff(adjsurv, loess_smoother=TRUE)

# plot without the line at reference
plot_curve_diff(adjsurv, line_at_ref=FALSE)
plot_curve_ratio(adjsurv, line_at_ref=FALSE)

# plot with some custom parameters
plot_curve_diff(adjsurv, conf_int=TRUE, color="blue", linetype="dotted",
                alpha=0.8, line_at_ref_size=1.1, line_at_ref_color="red",
                loess_smoother=TRUE, loess_span=0.55)

# adding a p-value for a difference test in the interval [0, 0.75]
plot_curve_diff(adjsurv, conf_int=TRUE, p_value=TRUE, integral_from=0,
                integral_to=0.75, integral=TRUE)

# adding a p-value for a difference test in the interval [0, 0.75],
# and also showing that integral visually in the plot
plot_curve_diff(adjsurv, conf_int=FALSE, p_value=TRUE, integral_from=0,
                integral_to=0.75, integral=TRUE, fill_area=TRUE,
                interval=TRUE)
}
}
