# Adjusted Survival Curves for Categorical Confounders using the Method by Amato (1988)

This page explains the details of estimating confounder-adjusted
survival curves using a weighted average of stratified Kaplan-Meier
estimates using the method described in Amato (1988)
(`method="strat_amato"` in the
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

- reference:

  A `data.frame` to be used as a reference population when weighting the
  survival curves or `NULL` (default). If `NULL` the survival curves are
  weighted in reference to the full sample supplied using `data`,
  regardless of the `variable` level. If a `data.frame` is supplied it
  needs to include all variables specified in `adjust_vars`.

## Details

- **Type of Adjustment:** The survival curves are adjusted by
  calculating a weighted version of the Kaplan-Meier estimator, based on
  stratification on covariates. This only works for categorical
  confounders. See below for more information.

- **Doubly-Robust:** Estimates are not Doubly-Robust.

- **Categorical groups:** Any number of levels in `variable` are
  allowed. Must be a factor variable.

- **Approximate Variance:** Calculations to approximate the variance and
  confidence intervals are not available. Bootstrap confidence intervals
  can however be calculated with all supported models. See
  [`?adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md)
  for more information on bootstrapping.

- **Allowed Time Values:** Allows both continuous and integer time.

- **Bounded Estimates:** Estimates are guaranteed to be bounded in the 0
  to 1 probability range.

- **Monotone Function:** Estimates are guaranteed to be monotone.

- **Dependencies:** This method has no dependencies.

This is one of the older adjustment methods described in the literature.
It only works for categorical confounders. If adjustments for continuous
confounders are desired, the user needs to explicitly categorize the
continuous confounders. It is recommended to use one of the other
methods implemented in this package in that case. The method works
exactly as described in Amato (1988). The number of people at risk and
the number of events in each stratum at each point in time is reweighted
and combined into a single estimate for each treatment. The reference
data used to calculate the weights is the pooled sample (`data`) by
default, but external reference data can be supplied. A more detailed
description can be found in the original article.

If a character vector is supplied in the `adjust_vars` argument, every
possible combination of the variables specified in `adjust_vars` will be
used as strata. This might lead to problems when there are strata with
very little data in them. In contrast to other stratification based
methods however, this method allows the estimation of adjusted survival
curves up to the last point in time at which there is at least one
individual at risk in the pooled sample.

## Value

Adds the following additional objects to the output of the
`adjustedsurv` function:

- `Pjs`: The weights used for each stratum.

## References

David A. Amato (1988). "A Generalized Kaplan-Meier Estimator for
Heterogenous Populations". In: Communications in Statistics: Theory and
Methods 17.1, pp. 263-286

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
                        method="strat_amato",
                        adjust_vars=c("x1", "x3"),
                        conf_int=FALSE)

# plot the curves
plot(adjsurv)
#> Ignoring unknown labels:
#> • linetype : "Group"
#> • fill : "Group"
```
