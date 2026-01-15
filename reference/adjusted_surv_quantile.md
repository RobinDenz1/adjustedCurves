# Estimate Confounder-Adjusted Survival Time Quantiles

This function can be utilized to estimate confounder-adjusted survival
time quantiles, including the median survival time, given previously
estimated adjusted survival curves. May also be used to estimate the
difference between or ratio of two adjusted survival time quantiles.

## Usage

``` r
adjusted_surv_quantile(adjsurv, p=0.5, conf_int=FALSE,
                       conf_level=adjsurv$conf_level,
                       use_boot=FALSE, interpolation="steps",
                       contrast="none", group_1=NULL,
                       group_2=NULL)
```

## Arguments

- adjsurv:

  An `adjustedsurv` object created using the `adjustedsurv` function.

- p:

  The quantile of interest. To calculate the median survival time, set
  this parameter to 0.5 (default). Multiple values in form of a numeric
  vector are allowed.

- conf_int:

  Whether to calculate confidence intervals or not. Those are calculated
  in the same way as the quantiles but using the confidence limit
  curves. This requires either that `conf_int=TRUE` or `bootstrap=TRUE`
  was used in the original `adjustedsurv` function call.

- conf_level:

  A single number specifying the confidence level of the confidence
  intervals. By default the previously estimated confidence intervals
  and thereby the same confidence level that was used in the original
  `adjustedsurv` call will be used. If `conf_level` is set to a
  different value, the confidence intervals will be re-calculated first.
  This works for all methods except `method="km"`.

- use_boot:

  Whether to use the bootstrap confidence interval estimates of the
  survival curves to estimate the confidence intervals of the survival
  time quantiles or not. Can only be used when `bootstrap=TRUE` was used
  in the original `adjustedsurv` function call. Ignored if
  `conf_int=FALSE`.

- interpolation:

  Either `"steps"` (default) or `"linear"`. This parameter controls how
  interpolation is performed. If this argument is set to `"steps"`, the
  curves will be treated as step functions. If it is set to `"linear"`,
  the curves wil be treated as if there are straight lines between the
  point estimates instead. Points that lie between estimated points will
  be interpolated accordingly.

- contrast:

  A single character string, specifying which contrast should be
  estimated. Needs to be one of `"none"` (estimate no contrasts, just
  return the adjusted survival time quantile, the default), `"diff"`
  (estimate the difference) or `"ratio"` (estimate the ratio). When
  `conf_int=TRUE` is also specified and bootstrapping was performed in
  the original `adjustedsurv` call, this function will also estimate the
  corresponding standard error, the confidence interval and a p-value
  testing whether the difference is equal to 0 (or the ratio is equal to
  1). To specify which difference/ratio should be calculated, the
  `group_1` and `group_2` arguments can be used. By default, the
  difference/ratio between the first and second level in `variable` is
  computed.

- group_1:

  Optional argument to get a specific difference or ratio. This argument
  takes a single character string specifying one of the levels of the
  `variable` used in the original `adjustedsurv` function call. This
  group will be subtracted from. For example if `group_1="A"` and
  `group_2="B"` and `contrast=="diff"` the difference `A - B` will be
  used. If `NULL`, the order of the factor levels in the original data
  determines the order. Ignored if `contrast="none"`.

- group_2:

  Also a single character string specifying one of the levels of
  `variable`. This corresponds to the right side of the difference/ratio
  equation. See argument `group_2`. Ignored if `contrast="none"`.

## Details

The median survival time is simply the time when half the patients are
expected to be alive. The chance of surviving beyond that time is 50
percent. In general, any quantile can be calculated this way. Those can
be read directly from the respective survival curve by drawing a
straight line from the desired quantile `p` on the Y-Axis and reading
the X-Axis value where this line intersects with the survival curve. The
adjusted survival time quantile for group \\z\\ (corresponding to the
`variable` argument in the
[`adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md)
function) is formally defined as:

\$\$\hat{Q}\_z(p) = min\left(t \| \hat{S}\_z(t) \leq p\right)\$\$

where \\\hat{S}\_z(t)\\ is the estimated counterfactual survival
function for \\z\\.

If the survival probability never drops below `p`, the survival time
quantile cannot be calculated. This function calculates this quantity
automatically.

***Confidence intervals***

When estimating the adjusted survival time quantiles directly, the
confidence intervals are simply read off the survival curves similarly
to the point estimates. For example, this means that the lower limit of
the confidence interval for some survival time quantile is simply the
first point in time at which the curve defined by the lower confidence
limit goes below \\p\\.

If `contrast="none"`, a different approach is used. Note that this
functionality works only with bootstrapping (regardless of whether
`use_boot` is `TRUE` or `FALSE`). In this case, the survival time
quantiles are calculated for each bootstrap sample separately first.
Afterwards, the standard deviation of the resulting estimates is
calculated for each group. For differences, the pooled standard error is
then estimated as \\SE\_{group_1 - group_2} = \sqrt{SE\_{group_1}^2 +
SE\_{group_2}^2}\\ and used along with the normal approximation to
calculate confidence intervals. This is similar to the strategy of Wu
(2011). The p-value is also estimated using the pooled standard error
and a one-sample two-sided t-test with the null-hypothesis that the
difference is equal to 0. The confidence interval and p-value of the
ratio is estimated using the method of Fieller (1954).

If one or both of the survival time quantiles used for the difference /
ratio calculation could not be estimated from the respective bootstrap
sample (due to the curves not extending that far) it is simply
discarded.

***More than Two Groups***

If more than two groups are present in `variable`, all other comparisons
except for `group_1 vs. group_2` are ignored. If multiple comparisons
are desired, the user needs to call this function multiple times and
adjust the `group_1` and `group_2` arguments accordingly.

***Multiple Imputation***

If multiple imputation was used in the original function call, the
survival time quantiles are read off the final pooled survival curves
directly. This also works when using `contrast="diff"` or
`contrast="ratio"`. When using either `contrast="diff"` or
`contrast="ratio"` in conjunction with `conf_int=TRUE` (plus
bootstrapping), the standard errors of the difference or ratio are
estimated for each imputed dataset separately and pooled using Rubins
Rule afterwards. The pooled standard errors are then used to perform the
same calculations as described above.

## Value

Returns a `data.frame` containing the columns `p` (the quantiles from
the original function call), `group` (groups in `variable`) and `q_surv`
(the survival time quantile).

If `conf_int=TRUE` was used it also includes the confidence limits in
the `ci_lower` and `ci_upper` columns.

If `contrast="diff"` was used, the resulting `data.frame` instead
contains the columns `p` (the quantiles from the original function all)
and `diff` (the difference between the two survival time quantiles). If
`conf_int=TRUE` was used it additionally contains the columns `se` (the
standard deviation of the bootstrap estimates), `ci_lower` and
`ci_upper` (the confidence interval) and `n_boot` (the number of
bootstrap samples that could be used). The output is similar when using
`contrast="ratio"`, except that the `diff` column is called `ratio`.

## References

Omer Ben-Aharon, Racheli Magnezi, Moshe Leshno, and Daniel A. Goldstein
(2019). "Median Survival or Mean Survival: Which Measure is the Most
Appropriate for Patients, Physicians, and Policymakers?" In: The
Oncologist 24, pp. 1469-1478

Zhongxue Chen and Guoyi Zhang (2016). "Comparing Survival Curves based
on Medians". In: BMC Medical Research Methodology 16.33

Jianrong Wu (2011). "Confidence Intervals for the Difference of Median
Failure Times Applied to Censored Tumor Growth Delay Data". In:
Statistics in Biopharmaceutical Research 3.3, pp. 488-496

Edgar C. Fieller (1954). "Some Problems in Interval Estimation". In:
Journal of the Royal Statistical Society, Series B 16.2, pp. 175-185

## Author

Robin Denz

## See also

[`adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md)

## Examples

``` r
library(adjustedCurves)
library(survival)

# simulate some data as example
set.seed(42)
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
                        bootstrap=FALSE)

# calculate adjusted median survival times
adjusted_surv_quantile(adjsurv)
#>     p group    q_surv
#> 1 0.5     0 0.5105218
#> 2 0.5     1 0.6561380

# calculate difference of adjusted median survival times between groups
adjusted_surv_quantile(adjsurv, contrast="diff")
#>     p       diff
#> 1 0.5 -0.1456162

# calculate other quantiles + confidence intervals
adjusted_surv_quantile(adjsurv, conf_int=TRUE, p=c(0.2, 0.4))
#>     p group    q_surv  ci_lower  ci_upper
#> 1 0.2     0 0.7408272 0.6177998 1.0331462
#> 2 0.4     0 0.6083660 0.4784817 0.7015162
#> 3 0.2     1        NA 0.7408272        NA
#> 4 0.4     1 0.7501622 0.6177998 1.0503161
```
