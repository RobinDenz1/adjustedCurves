# Estimate Confounder-Adjusted Restricted Mean Survival Times

This function can be utilized to estimate the confounder-adjusted
restricted mean survival time, given previously estimated adjusted
survival curves.

## Usage

``` r
adjusted_rmst(adjsurv, to, from=0, conf_int=FALSE,
              conf_level=0.95, interpolation="steps",
              contrast="none", group_1=NULL,
              group_2=NULL)
```

## Arguments

- adjsurv:

  An `adjustedsurv` object created using the `adjustedsurv` function.

- from:

  A single number specifying the left side of the time interval of
  interest. See details. Usually this should be kept at 0 (default) to
  estimate the standard RMST. Should only be changed if there are good
  reasons for it.

- to:

  One or multiple numbers specifying the right side of the time interval
  of interest. If a vector of numbers is supplied, the adjusted RMST
  will be estimated for each value of `to`. See details.

- conf_int:

  Whether bootstrap estimates should be used to estimate the standard
  errors and confidence intervals of the RMST estimates. Can only be
  used if `bootstrap=TRUE` was used in the
  [`adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md)
  call.

- conf_level:

  A number specifying the confidence level of the bootstrap confidence
  intervals.

- interpolation:

  Either `"steps"` (default) or `"linear"`. This parameter controls how
  interpolation is performed. If this argument is set to `"steps"`, the
  curves will be treated as step functions. If it is set to `"linear"`,
  the curves wil be treated as if there are straight lines between the
  point estimates instead. Points that lie between estimated points will
  be interpolated accordingly. Should usually be kept at `"steps"`. See
  Details.

- contrast:

  A single character string, specifying which contrast should be
  estimated. Needs to be one of `"none"` (estimate no contrasts, just
  return the adjusted RMST, the default), `"diff"` (estimate the
  difference) or `"ratio"` (estimate the ratio). When `conf_int=TRUE` is
  also specified and bootstrapping was performed in the original
  `adjustedsurv` call, this function will also estimate the
  corresponding standard error, the confidence interval and a p-value
  testing whether the difference is equal to 0 (or the ratio is equal to
  1). To specify which difference/ratio should be calculated, the
  `group_1` and `group_2` arguments can be used. By default, the
  difference/ratio between the first and second level in `variable` is
  computed.

- group_1:

  Optional argument to get a specific difference or ratio. This argument
  takes a single character string specifying one of the levels of the
  `variable` used in the original `adjustedsurv` or `adjustedcif`
  function call. This group will be subtracted from. For example if
  `group_1="A"` and `group_2="B"` and `contrast="diff"` the difference
  `A - B` will be used. If `NULL`, the order of the factor levels in the
  original `data` determines the order. Ignored if `contrast="none"`.

- group_2:

  Also a single character string specifying one of the levels of
  `variable`. This corresponds to the right side of the difference/ratio
  equation. See argument `group_2`. Ignored if `contrast="none"`.

## Details

The adjusted restricted mean survival times (RMST) are estimated by
integrating the estimated adjusted survival curves in a specified
interval. Let \\Z\\ be the grouping variable (corresponding to the
`variable` argument in the
[`adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md)
function) with possible levels \\Z \in \\0, 1, 2, ..., k\\\\. \\T\\ is
defined as the time and \\\hat{S}\_z(t)\\ denotes the estimated
counterfactual survival function. The RMST is then defined as:

\$\$RMST\_{z}(to) = \int\_{from=0}^{to} \hat{S}\_z(t)dt\$\$

It can be interpreted as the mean survival time of individuals in group
\\Z = z\\ in the interval \[`from`, `to`\]. Note however that simply
subtracting the estimates from each other does not give a correct
estimate of the area between the survival curves if the respective
curves cross at some point. The
[`adjusted_curve_test`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_curve_test.md)
function can be used to calculate the actual area between the curves
instead. See
[`?adjusted_curve_test`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_curve_test.md)
for more information.

***Confidence Intervals***

If the `adjsurv` object was created with `bootstrap=TRUE` in the
[`adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md)
function, bootstrap confidence intervals and standard errors for the
RMSTs can be approximated by setting `conf_int` to `TRUE`. If bootstrap
samples occur where the survival function is not estimated up to `to`,
the bootstrap sample is discarded and not used in further calculations.
Approximate variance calculations not relying on the bootstrap estimates
are currently not implemented. When using `contrast="diff"` the standard
error of the difference between the two RMST values is approximated by
\\SE\_{group_1 - group_2} = \sqrt{SE\_{group_1}^2 + SE\_{group_2}^2}\\.
When using `contrast="ratio"` the confidence intervals are calculated
using the approximate formula given by Fieller (1954), assuming that the
values are independent.

***More than Two Groups***

If more than two groups are present in `variable`, all other comparisons
except for `group_1 vs. group_2` are ignored. If multiple comparisons
are desired, the user needs to call this function multiple times and
adjust the `group_1` and `group_2` arguments accordingly.

***Multiple Imputation***

If multiple imputation was used when creating the `adjsurv` object, the
analysis is carried out on all multiply imputed datasets and pooled
using Rubins Rule. When bootstrapping was carried out as well, the
pooled standard error over all imputed datasets is used in combination
with the normal approximation to re-calculate the bootstrap confidence
intervals.

***Competing Risks***

This function cannot be used with `adjustedcif` objects, because the
survival probability cannot be estimated in an unbiased way when
competing risks are present. However, a very similar quantity, the
*adjusted restricted mean time lost*, can be calculated using the
[`adjusted_rmtl`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_rmtl.md)
function.

***Graphical Displays***

A plot of the RMST over time (with changing values for the `to`
argument) can be produced using the
[`plot_rmst_curve`](https://robindenz1.github.io/adjustedCurves/reference/plot_rmst_curve.md)
function.

***Computational Details***

Instead of relying on numerical integration, this function uses exact
calculations. This is achieved by using either step-function
interpolation (`interpolation="steps"`, the default) or linear
interpolation (`interpolation="linear"`). In the former case, the
integral is simply the sum of the area of the squares defined by the
step function. In the second case, the integral is simply the sum of the
area of the rectangles. Either way, there is no need for approximations.
In some situations (for example when using parametric survival models
with `method="direct"`), the curves are not step functions. In this case
the `interpolation` argument should be set to `"linear"`.

## Value

Returns a `data.frame` containing the columns `to` (the supplied `to`
values), `group` (groups in `variable`) and `rmst` (the estimated
restricted mean survival time).

If `conf_int=TRUE` was used it additionally contains the columns `se`
(the standard error of the restricted mean survival time), `ci_lower`
(lower limit of the confidence interval), `ci_upper` (upper limit of the
confidence interval) and `n_boot` (the actual number of bootstrap
estimates used).

If `contrast="diff"` was used, it instead returns a `data.frame` that
contains the columns `to`, `diff` (the difference between the RMST
values), `se` (the standard error of the difference), `ci_lower` (lower
limit of the confidence interval), `ci_upper` (upper limit of the
confidence interval) and `p_value` (the p-value of the one-sample
t-test). The same results are presented when using `contrast="ratio"`,
except that the `diff` column is named `ratio` and that there is no `se`
column.

## References

Sarah C. Conner, Lisa M. Sullivan, Emelia J. Benjamin, Michael P.
LaValley, Sandro Galea, and Ludovic Trinquart (2019). "Adjusted
Restricted Mean Survival Times in Observational Studies". In: Statistics
in Medicine 38, pp. 3832-3860

Patrick Royston and Mahesh K. B. Parmar (2013). "Restricted Mean
Survival Time: An Alternative to the Hazard Ratio for the Design and
Analysis of Randomized Trials with a Time-To-Event Outcome". In: BMC
Medical Research Methodology 13.152

Edgar C. Fieller (1954). "Some Problems in Interval Estimation". In:
Journal of the Royal Statistical Society, Series B 16.2, pp. 175-185

## Author

Robin Denz

## See also

[`adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md),
[`plot_rmst_curve`](https://robindenz1.github.io/adjustedCurves/reference/plot_rmst_curve.md)

## Examples

``` r
library(adjustedCurves)
library(survival)

set.seed(42)

# simulate some data as example
sim_dat <- sim_confounded_surv(n=30, max_t=1.2)
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
                        conf_int=FALSE,
                        bootstrap=TRUE,
                        n_boot=10) # n_boot should be much higher in reality

# calculate adjusted restricted mean survival times from 0 to 0.5
adjrmst <- adjusted_rmst(adjsurv, from=0, to=0.5, conf_int=FALSE)

# calculate adjusted restricted mean survival times from 0 to 0.5
# and from 0 to 0.7 simultaneously
adjrmst <- adjusted_rmst(adjsurv, from=0, to=c(0.5, 0.7), conf_int=FALSE)

# calculate adjusted restricted mean survival times from 0 to 0.5,
# including standard errors and confidence intervals
adjrmst <- adjusted_rmst(adjsurv, from=0, to=0.5, conf_int=TRUE,
                         conf_level=0.95)

# calculate difference between adjusted restricted mean survival times
# from 0 to 0.5 in the two groups
adjrmst <- adjusted_rmst(adjsurv, from=0, to=0.5, conf_int=FALSE,
                         contrast="diff")

# calculate ratio between adjusted restricted mean survival times
# from 0 to 0.5 in the two groups
adjrmst <- adjusted_rmst(adjsurv, from=0, to=0.5, conf_int=FALSE,
                         contrast="ratio")
```
