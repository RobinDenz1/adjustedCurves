# Inverse Probability of Treatment Weighted Kaplan-Meier estimates

This page explains the details of estimating inverse probability of
treatment weighted survival curves using a weighted version of the
Kaplan-Meier estimator for single event time-to-event data
(`method="iptw_km"` in the
[`adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md)
function). All regular arguments of the `adjustedsurv` function can be
used. Additionally, the `treatment_model` argument has to be specified
in the `adjustedsurv` call. Further arguments specific to this method
are listed below.

## Arguments

- treatment_model:

  \[**required**\] Must be either a model object with `variable` as
  response variable, a vector of weights or a formula which can be
  passed to `WeightIt`.

- weight_method:

  Method used in `WeightIt` function call. Ignored if `treatment_model`
  is not a formula object. Defaults to `"ps"`.

- stabilize:

  Whether to stabilize the weights or not. Is set to `FALSE` by default.
  Stabilizing weights ensures that the sum of all weights is equal to
  the original sample size. It has no effect on point estimates, only on
  the asymptotic variance calculations and confidence intervals.

- trim:

  Can be either `FALSE` (default) or a numeric value at which to trim
  the weights. If `FALSE`, weights are used as calculated or supplied.
  If a numeric value is supplied, all weights that are bigger than
  `trim` are set to `trim` before the analysis is carried out. Useful
  when some weights are extremely large.

- trim_quantiles:

  Alternative argument to trim weights based on quantiles. Can be either
  `FALSE` (default) to use no trimming, or a numeric vector containing
  exactly two values between 0 and 1. These values specify the quantiles
  that the weights should be trimmed at. For example, if `c(0.01, 0.99)`
  is supplied to this argument, all weights that are lower than the 0.01
  quantile of the weight distribution will be set to that quantile and
  all weights that are higher than the 0.99 quantile of the weight
  distributions will be set to the 0.99 quantile.

- extend_to_last:

  Either `TRUE` (default) or `FALSE`, indicating whether the survival
  curves should be extended up to the last censored observation time
  (`TRUE`) or only to the last event time (`FALSE`).

- ...:

  Further arguments passed to
  [`weightit`](https://ngreifer.github.io/WeightIt/reference/weightit.html).

## Details

- **Type of Adjustment:** Requires a model describing the treatment
  assignment mechanism. This must be either a
  [`glm`](https://rdrr.io/r/stats/glm.html) or
  [`multinom`](https://rdrr.io/pkg/nnet/man/multinom.html) object.

- **Doubly-Robust:** Estimates are not Doubly-Robust.

- **Categorical groups:** Any number of levels in `variable` are
  allowed. Must be a factor variable.

- **Approximate Variance:** Calculations to approximate the variance and
  confidence intervals are available.

- **Allowed Time Values:** Allows both continuous and integer time.

- **Bounded Estimates:** Estimates are guaranteed to be bounded in the 0
  to 1 probability range.

- **Monotone Function:** Estimates are guaranteed to be monotone.

- **Dependencies:** This method does not depend on other packages
  directly. However the WeightIt package is required if
  `treatment_model` is a formula object.

This method works by modeling the treatment assignment mechanism.
Adjusted survival curves are calculated by first estimating appropriate
case-weights for each observation in `data`. This can be done using
inverse probability of treatment weights using the propensity score
(usually estimated using a logistic regression model) or by some other
method (see
[`?weightit`](https://ngreifer.github.io/WeightIt/reference/weightit.html)).
Those weights are used in a weighted version of the Kaplan-Meier
estimator proposed by Xie and Liu (2005). If the weights are correctly
estimated the resulting estimates will be unbiased. The only difference
to the
[`iptw_cox`](https://robindenz1.github.io/adjustedCurves/reference/surv_iptw_cox.md)
method is a slightly different weighting approach.

Asymptotic variances are calculated using the equations given in Xie and
Liu (2005). It is also recommended to use stabilized weights by using
`stabilize=TRUE` (the default value). More information can be found in
the cited literature.

## Value

Adds the following additional objects to the output of the
`adjustedsurv` function:

- `weights`: The final weights used in the analysis.

- `n_at_risk`: A `data.frame` containing the weighted number at risk and
  weighted number of events used in the calculations at each point in
  time for both groups.

## References

Jun Xie and Chaofeng Liu (2005). "Adjusted Kaplan-Meier Estimator and
Log- Rank Test with Inverse Probability of Treatment Weighting for
Survival Data". In: Statistics in Medicine 24, pp. 3089-3110

Stanley Xu, Colleen Ross and Marsha A. Raebel, Susan Shetterly,
Christopher Blanchette, and David Smith (2010). "Use of Stabilized
Inverse Propensity Scores as Weights to Directly Estimate Relative Risk
and Its Confidence Intervals". In: Value in Health 13.2, pp. 273-277

## Author

Robin Denz

## See also

[`weightit`](https://ngreifer.github.io/WeightIt/reference/weightit.html)

## Examples

``` r
library(adjustedCurves)

set.seed(42)

# simulate some data as example
sim_dat <- sim_confounded_surv(n=50, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# estimate a treatment assignment model
glm_mod <- glm(group ~ x1 + x3 + x5 + x6, data=sim_dat, family="binomial")

# use it to calculate adjusted survival curves
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="iptw_km",
                        treatment_model=glm_mod)

# Alternatively, use custom weights
# In this example we use weights calculated using the propensity score,
# which is equal to using the glm model directly in the function
ps_score <- glm_mod$fitted.values
weights <- ifelse(sim_dat$group==1, 1/ps_score, 1/(1-ps_score))

adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="iptw_km",
                        treatment_model=weights)

if (requireNamespace("WeightIt")) {

# And a third alternative: use the WeightIt package
# here an example with equal results to the ones above:
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="iptw_km",
                        treatment_model=group ~ x1 + x3 + x5 + x6,
                        weight_method="ps")

# here an example using Entropy Balancing Weighting:
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="iptw_km",
                        treatment_model=group ~ x1 + x3 + x5 + x6,
                        weight_method="ebal")
}
```
