\name{plot_rmst_curve}
\alias{plot_rmst_curve}

\title{
Plot Adjusted Restricted Mean Survival Time Curves
}
\description{
A function to graphically display the Restricted Mean Survival Time (RMST) over time, using confounder-adjusted survival curves which where previously estimated using the \code{\link[adjustedCurves]{adjustedsurv}} function. As the other plot functions in this package, it internally uses the \pkg{ggplot2} package and allows a variety of options. Alternatively plots the difference or ratio between two RMST Curves.
}
\usage{
plot_rmst_curve(adjsurv, times=NULL, conf_int=FALSE,
                conf_level=0.95, interpolation="steps",
                contrast="none", group_1=NULL, group_2=NULL,
                max_t=Inf, color=TRUE, linetype=FALSE,
                facet=FALSE, size=1, alpha=1, xlab="Time",
                ylab="RMST", title=NULL, subtitle=NULL,
                legend.title="Group", legend.position="right",
                gg_theme=ggplot2::theme_classic(),
                custom_colors=NULL, custom_linetypes=NULL,
                conf_int_alpha=0.4,
                line_at_ref=TRUE, line_at_ref_size=0.7,
                line_at_ref_color="grey",
                line_at_ref_linetype="dashed",
                line_at_ref_alpha=1, ...)
}

\arguments{
  \item{adjsurv}{
An \code{adjustedsurv} object created using the \code{\link[adjustedCurves]{adjustedsurv}} function.
  }
  \item{times}{
A vector of points in time, passed to the \code{to} argument of the \code{\link{adjusted_rmst}} function or \code{NULL} (default). If \code{NULL}, the adjusted RMST is estimated at all points at which an event occurred. Otherwise it is estimated at \code{times}.
  }
  \item{conf_int}{
A logical variable indicating whether the bootstrap confidence intervals should be drawn.
  }
  \item{conf_level}{
Corresponds to the argument of the same name in the \code{\link{adjusted_rmst}} function.
  }
  \item{interpolation}{
Corresponds to the argument of the same name in the \code{\link{adjusted_rmst}} function.
  }
  \item{contrast}{
Which contrast between two adjusted RMST curves should be plotted. Should be one of \code{c("none", "diff", "ratio")}. See argument of the same name in the \code{\link{adjusted_rmst}} function.
  }
  \item{group_1}{
A single character string specifying one of the possible levels of \code{variable}. This can be used in conjunction with the \code{group_2} argument to control how the difference or ratio should be calculated when using either \code{contrast="diff"} or \code{contrast="ratio"}. Ignored when \code{contrast="none"} (default). See argument of the same name in the \code{\link{adjusted_rmst}} function.
  }
  \item{group_2}{
See \code{group_2}.
  }
  \item{max_t}{
A number indicating the latest survival time which is to be plotted.
  }
  \item{color}{
A logical variable indicating whether the curves should be colored differently. The \code{custom_colors} argument can be used to directly specify which colors to use. Set to \code{FALSE} to keep the plot black and white. If \code{contrast="diff"} or \code{contrast="ratio"} are used, a character string specifying a color that should be used for the plot can be supplied to this argument directly.
  }
  \item{linetype}{
A logical variable indicating whether the curves should have different linetypes. The \code{custom_linetypes} argument can be used to directly specify which linetypes to use. Set to \code{FALSE} to keep all lines solid. If \code{contrast="diff"} or \code{contrast="ratio"} are used, a character string specifying a linetype that should be used for the plot can be supplied to this argument directly.
  }
  \item{facet}{
A logical variable indicating whether the curves should be in different facets.
  }
  \item{size}{
A number controlling the thickness of the RMST curves.
  }
  \item{alpha}{
A number controlling the transparency level of the RMST curves.
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
  \item{subtitle}{
A character string to be used as the subtitle of the plot. Set to \code{NULL} if no subtitle should be used.
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
  \item{custom_colors}{
A (named) vector to specify the colors of each adjusted RMST curve and possibly its confidence region. Set to \code{NULL} to use the \code{ggplot2} default values. Ignored if \code{color=FALSE}.
  }
  \item{custom_linetypes}{
A (named) vector to specify the linetype of each adjusted RMST curve. Set to \code{NULL} to use the \code{ggplot2} default values. Ignored if \code{color=FALSE}. Ignored if \code{linetype=FALSE}.
  }
  \item{conf_int_alpha}{
A number indicating the level of transparency that should be used when drawing the confidence regions.
  }
  \item{line_at_ref}{
Whether to draw a line at the reference value. This line is drawn at 0 if \code{contrast="diff"} and at 1 if \code{contrast="ratio"}. This and all associated argument are ignored otherwise.
  }
  \item{line_at_ref_size}{
A single number specifying the thickness of the line at the reference value.
  }
  \item{line_at_ref_color}{
A single character string specifying the color of the line at the reference value.
  }
  \item{line_at_ref_linetype}{
A single character string specifying the linetype of the line at the reference value.
  }
  \item{line_at_ref_alpha}{
A single number between 0 and 1 specifying the transparency level of the line at the reference value.
  }
  \item{...}{
Currently not used.
  }
}
\details{
This function simply calls the \code{\link{adjusted_rmst}} for a range of \code{to} values, getting adjusted RMST estimates over the whole range of the survival curves. Those estimates are then plotted as a curve with the adjusted RMST replacing the survival probability on the Y-Axis. For a brief description on the RMST and how it is calculated in this package, see the documentation of the \code{\link{adjusted_rmst}} function. Literature describing the RMST Curve Plots in more detail is given in the references section.

The RMST curve can only be created for adjusted survival curves. A similar graphic for the adjusted CIFs can be created by utilizing the adjusted Restricted Mean Time Lost (RMTL). The calculation of that statistic is implemented in the \code{\link{adjusted_rmtl}} function and the associated curve can be created using the \code{\link{plot_rmtl_curve}} function.

If confidence intervals are specified and there are many points in time in \code{times}, this function might get somewhat slow. It will be even slower if multiple imputation was also used when creating the \code{adjustedsurv} object.
}
\value{
Returns a \code{ggplot2} object.
}
\references{
Lihui Zhao, Brian Claggett, Lu Tian, Hajime Uno, Marc A. Pfeffer, Scott D. Solomon, Lorenzo Trippa, and L. J. Wei (2016). "On the Restricted Mean Survival Time Curve in Survival Analysis". In: Biometrics 72.1, pp. 215-221

Jason J. Z. Liao, Frank Liu, and Wen-Chi Wu (2020). "Dynamic RMST Curves for Survival Analysis in Clinical Trials". In: BMC Medical Research Methodology 20.218
}
\author{
Robin Denz
}
\seealso{
\code{\link{adjustedsurv}}, \code{\link{adjusted_rmst}}, \code{\link[ggplot2]{ggplot}}
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

# plot the curves with default values
plot_rmst_curve(adjsurv)

# plot with confidence intervals
plot_rmst_curve(adjsurv, conf_int=TRUE)

# plot the difference instead
plot_rmst_curve(adjsurv, contrast="diff")

# plot with some custom options
plot_rmst_curve(adjsurv, max_t=0.5, linetype=TRUE,
                custom_colors=c("green", "blue"))
}
}
