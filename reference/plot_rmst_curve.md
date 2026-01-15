# Plot Adjusted Restricted Mean Survival Time Curves

A function to graphically display the Restricted Mean Survival Time
(RMST) over time, using confounder-adjusted survival curves which where
previously estimated using the
[`adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md)
function. As the other plot functions in this package, it internally
uses the ggplot2 package and allows a variety of options. Alternatively
plots the difference or ratio between two RMST Curves.

## Usage

``` r
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
```

## Arguments

- adjsurv:

  An `adjustedsurv` object created using the
  [`adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md)
  function.

- times:

  A vector of points in time, passed to the `to` argument of the
  [`adjusted_rmst`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_rmst.md)
  function or `NULL` (default). If `NULL`, the adjusted RMST is
  estimated at all points at which an event occurred. Otherwise it is
  estimated at `times`.

- conf_int:

  A logical variable indicating whether the bootstrap confidence
  intervals should be drawn.

- conf_level:

  Corresponds to the argument of the same name in the
  [`adjusted_rmst`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_rmst.md)
  function.

- interpolation:

  Corresponds to the argument of the same name in the
  [`adjusted_rmst`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_rmst.md)
  function.

- contrast:

  Which contrast between two adjusted RMST curves should be plotted.
  Should be one of `c("none", "diff", "ratio")`. See argument of the
  same name in the
  [`adjusted_rmst`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_rmst.md)
  function.

- group_1:

  A single character string specifying one of the possible levels of
  `variable`. This can be used in conjunction with the `group_2`
  argument to control how the difference or ratio should be calculated
  when using either `contrast="diff"` or `contrast="ratio"`. Ignored
  when `contrast="none"` (default). See argument of the same name in the
  [`adjusted_rmst`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_rmst.md)
  function.

- group_2:

  See `group_2`.

- max_t:

  A number indicating the latest survival time which is to be plotted.

- color:

  A logical variable indicating whether the curves should be colored
  differently. The `custom_colors` argument can be used to directly
  specify which colors to use. Set to `FALSE` to keep the plot black and
  white. If `contrast="diff"` or `contrast="ratio"` are used, a
  character string specifying a color that should be used for the plot
  can be supplied to this argument directly.

- linetype:

  A logical variable indicating whether the curves should have different
  linetypes. The `custom_linetypes` argument can be used to directly
  specify which linetypes to use. Set to `FALSE` to keep all lines
  solid. If `contrast="diff"` or `contrast="ratio"` are used, a
  character string specifying a linetype that should be used for the
  plot can be supplied to this argument directly.

- facet:

  A logical variable indicating whether the curves should be in
  different facets.

- size:

  A number controlling the thickness of the RMST curves.

- alpha:

  A number controlling the transparency level of the RMST curves.

- xlab:

  A character string to be used as the X-Axis label of the plot.

- ylab:

  A character string to be used as the Y-Axis label of the plot.

- title:

  A character string to be used as the title of the plot. Set to `NULL`
  if no title should be used.

- subtitle:

  A character string to be used as the subtitle of the plot. Set to
  `NULL` if no subtitle should be used.

- legend.title:

  A character string to be used as the title of the legend. Set to
  `NULL` if no legend should be included.

- legend.position:

  A character string specifying the position of the legend. Ignored if
  `legend_title=NULL`.

- gg_theme:

  A `ggplot2` theme object which will be used for the plot.

- custom_colors:

  A (named) vector to specify the colors of each adjusted RMST curve and
  possibly its confidence region. Set to `NULL` to use the `ggplot2`
  default values. Ignored if `color=FALSE`.

- custom_linetypes:

  A (named) vector to specify the linetype of each adjusted RMST curve.
  Set to `NULL` to use the `ggplot2` default values. Ignored if
  `color=FALSE`. Ignored if `linetype=FALSE`.

- conf_int_alpha:

  A number indicating the level of transparency that should be used when
  drawing the confidence regions.

- line_at_ref:

  Whether to draw a line at the reference value. This line is drawn at 0
  if `contrast="diff"` and at 1 if `contrast="ratio"`. This and all
  associated argument are ignored otherwise.

- line_at_ref_size:

  A single number specifying the thickness of the line at the reference
  value.

- line_at_ref_color:

  A single character string specifying the color of the line at the
  reference value.

- line_at_ref_linetype:

  A single character string specifying the linetype of the line at the
  reference value.

- line_at_ref_alpha:

  A single number between 0 and 1 specifying the transparency level of
  the line at the reference value.

- ...:

  Currently not used.

## Details

This function simply calls the
[`adjusted_rmst`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_rmst.md)
for a range of `to` values, getting adjusted RMST estimates over the
whole range of the survival curves. Those estimates are then plotted as
a curve with the adjusted RMST replacing the survival probability on the
Y-Axis. For a brief description on the RMST and how it is calculated in
this package, see the documentation of the
[`adjusted_rmst`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_rmst.md)
function. Literature describing the RMST Curve Plots in more detail is
given in the references section.

The RMST curve can only be created for adjusted survival curves. A
similar graphic for the adjusted CIFs can be created by utilizing the
adjusted Restricted Mean Time Lost (RMTL). The calculation of that
statistic is implemented in the
[`adjusted_rmtl`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_rmtl.md)
function and the associated curve can be created using the
[`plot_rmtl_curve`](https://robindenz1.github.io/adjustedCurves/reference/plot_rmtl_curve.md)
function.

If confidence intervals are specified and there are many points in time
in `times`, this function might get somewhat slow. It will be even
slower if multiple imputation was also used when creating the
`adjustedsurv` object.

## Value

Returns a `ggplot2` object.

## References

Lihui Zhao, Brian Claggett, Lu Tian, Hajime Uno, Marc A. Pfeffer, Scott
D. Solomon, Lorenzo Trippa, and L. J. Wei (2016). "On the Restricted
Mean Survival Time Curve in Survival Analysis". In: Biometrics 72.1, pp.
215-221

Jason J. Z. Liao, Frank Liu, and Wen-Chi Wu (2020). "Dynamic RMST Curves
for Survival Analysis in Clinical Trials". In: BMC Medical Research
Methodology 20.218

## Author

Robin Denz

## See also

[`adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md),
[`adjusted_rmst`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_rmst.md),
[`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)

## Examples

``` r
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
#> Ignoring unknown labels:
#> • fill : "Group"
```
