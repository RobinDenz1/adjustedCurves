# Adjusted Survival Curves for Categorical Confounders using the Method by Gregory (1988) and Nieto & Coresh (1996)

This page explains the details of estimating confounder-adjusted
survival curves using a weighted average of stratified Kaplan-Meier
estimates using the method described in Gregory (1988) and Nieto &
Coresh (1996) (`method="strat_gregory_nieto"` in the
[`adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md)
function). All regular arguments of the `adjustedsurv` function can be
used. Additionally, the `adjust_vars` argument has to be specified in
the `adjustedsurv` call. Further arguments specific to this method are
listed below.

## Arguments

- adjust_vars:

  \[**required**\] A single string or character vector specifying column
  names in data for which the survival curves should be adjusted for.
  The variables specified can be integers, factors or characters. Only
  categorical variables can be used with this method. See details.

## Details

- **Type of Adjustment:** The survival curves are adjusted by taking a
  weighted average of stratified Kaplan-Meier estimates. This only works
  for categorical confounders. See below for more information.

- **Doubly-Robust:** Estimates are not Doubly-Robust.

- **Categorical groups:** Any number of levels in `variable` are
  allowed. Must be a factor variable.

- **Approximate Variance:** Calculations to approximate the variance and
  confidence intervals are available. The estimator for the variance can
  be found in the appendix of Nieto & Coresh (1996).

- **Allowed Time Values:** Allows both continuous and integer time.

- **Bounded Estimates:** Estimates are guaranteed to be bounded in the 0
  to 1 probability range.

- **Monotone Function:** Estimates are guaranteed to be monotone.

- **Dependencies:** This method has no additional dependencies.

This is one of the older adjustment methods described in the literature.
It only works for categorical confounders. If adjustments for continuous
confounders are desired, the user needs to explicitly categorize the
continuous confounders. It is recommended to use one of the other
methods implemented in this package in that case. The method works
exactly as described in Gregory (1988). Similarly to the method
described in
[strat_cupples](https://robindenz1.github.io/adjustedCurves/reference/surv_strat_cupples.md),
Kaplan-Meier estimates are calculated for each strata and a weighted
average is taken. The only difference is a slightly different weighting
scheme. Weights are calculated using the pooled sample (`data`). In
contrast to other stratification based methods, external reference data
is not allowed. A more detailed description can be found in the original
article.

If a character vector is supplied in the `adjust_vars` argument, the
Kaplan-Meier estimates are created for each combination of all supplied
variables. If the sample size is small and/or there are many levels in
these variables, the estimates can become unstable or undefined. Because
it is a weighted average of Kaplan-Meier curves, estimates for this
method are only defined for points in time with a valid Kaplan-Meier
estimate in all strata. For example, if the Kaplan-Meier curve of the
strata "Treatment + male" only extends to t = 100, it will be impossible
to estimate the adjusted survival curve for t \> 100 using this method.

Nieto & Coresh (1996) proposed a very similar method. The only major
difference is that Nieto & Coresh (1996) used the control group as
reference population, which results in a different causal estimand.
Using the method by Nieto & Coresh (1996) with the full `data` as
reference population as described in Gregory (1988) produces exactly the
same results. Nieto & Coresh (1996) seemed to be unaware of the method
by Gregory (1988), as they did not mention it in their article. In
contrast to Gregory (1988) they however also proposed an approximate
estimator of the variance, which is implemented here. Their formulation
of this estimator also allows the use of time-dependent covariates and
left-truncated data. This is however not implemented here.

## Value

Adds no additional objects to the output of the `adjustedsurv` function.

## References

W. M. Gregory (1988). "Adjusting Survival Curves for Imbalances in
Prognostic Factors". In: British Journal of Cancer 58, pp. 202-204

F. Javier Nieto and Josef Coresh (1996). "Adjusting Survival Curves for
Confounders: A Review and a New Method". In: American Journal of
Epidemiology 143.10, pp. 1059-1068

## Author

Robin Denz

## See also

[`adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md)

## Examples

``` r
library(adjustedCurves)
library(survival)

set.seed(42)

# simulate some data as example
sim_dat <- sim_confounded_surv(n=50, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# adjust survival curves for some categorical confounders
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="strat_nieto",
                        adjust_vars=c("x1", "x3"),
                        conf_int=FALSE)

# plot the curves
plot(adjsurv)
#> Ignoring unknown labels:
#> • linetype : "Group"
#> • fill : "Group"
#> Warning: Removed 7 rows containing missing values or values outside the scale range
#> (`geom_step()`).
```
