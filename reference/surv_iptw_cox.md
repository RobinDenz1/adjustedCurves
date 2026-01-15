# Inverse Probability of Treatment Weighted Survival using Cox-Regression

This page explains the details of estimating inverse probability of
treatment weighted survival curves using a weighted univariate
cox-regression for single event time-to-event data (`method="iptw_cox"`
in the
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
  confidence intervals are not available.

- **Allowed Time Values:** Allows both continuous and integer time.

- **Bounded Estimates:** Estimates are guaranteed to be bounded in the 0
  to 1 probability range.

- **Monotone Function:** Estimates are guaranteed to be monotone.

- **Dependencies:** This method relies on the survival package.
  Additionally, the WeightIt package is required if `treatment_model` is
  a formula object.

This method works by modeling the treatment assignment mechanism.
Adjusted survival curves are calculated by first estimating appropriate
case-weights for each observation in `data`. This can be done using
inverse probability of treatment weights using the propensity score
(usually estimated using a logistic regression model) or by some other
method (see
[`?weightit`](https://ngreifer.github.io/WeightIt/reference/weightit.html)).
Those estimates are then used to fit a weighted Cox-Regression model,
stratified by `variable`. Survival Curves based on this model are
estimated using the method implemented in the `survfit.coxph` function.
More information can be found in the literature listed under
"references". The only difference to the
[`iptw_km`](https://robindenz1.github.io/adjustedCurves/reference/surv_iptw_km.md)
method is a slightly different weighting approach.

By default this method uses a a robust sandwich-type variance estimator
(`robust=TRUE` in the `coxph` function call) to calculate the standard
error used in the construction of confidence intervals. This estimator
has been shown to be biased by Austin (2016). Coupled with stabilized
weights however (`stabilize=TRUE`) this gives conservative estimates of
the variance and confidence intervals (Xu et al. 2010). It is still
recommended to use bootstrap confidence intervals instead. This can be
done by setting `bootstrap=TRUE` in the `adjustedsurv` function call.

## Value

Adds the following additional objects to the output of the
`adjustedsurv` function:

- `cox_model`: The stratified and weighted `coxph` model.

- `survfit`: The `survfit` object created using the `cox_model` object.

- `weights`: The final weights used in the analysis.

Returns a `list` object containing a `data.frame` with the estimated
adjusted survival probabilities for some points in time for each level
of `variable`, the weighted `coxph` model, the weighted `survfit` object
and the weights used in the analysis.

## References

Stephen R. Cole and Miguel A. Hernán (2004). "Adjusted Survival Curves
with Inverse Probability Weights". In: Computer Methods and Programs in
Biomedicine 2003.75, pp. 45-49

Peter C. Austin (2016). "Variance Estimation when Using Inverse
Probability of Treatment Weighting (IPTW) with Survival Analysis". In:
Statistics in Medicine 35, pp. 5642-5655

Stanley Xu, Colleen Ross and Marsha A. Raebel, Susan Shetterly,
Christopher Blanchette, and David Smith (2010). "Use of Stabilized
Inverse Propensity Scores as Weights to Directly Estimate Relative Risk
and Its Confidence Intervals". In: Value in Health 13.2, pp. 273-277

## Author

Robin Denz

## See also

[`weightit`](https://ngreifer.github.io/WeightIt/reference/weightit.html),
[`coxph`](https://rdrr.io/pkg/survival/man/coxph.html),
[`survfit.coxph`](https://rdrr.io/pkg/survival/man/survfit.coxph.html)

## Examples

``` r
library(adjustedCurves)
library(survival)

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
                        method="iptw_cox",
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
                        method="iptw_cox",
                        treatment_model=weights)

if (requireNamespace("WeightIt")) {

# And a third alternative: use the WeightIt package
# here an example with equal results to the ones above:
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="iptw_cox",
                        treatment_model=group ~ x1 + x3 + x5 + x6,
                        weight_method="ps")

# here an example using Entropy Balancing Weighting:
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="iptw_cox",
                        treatment_model=group ~ x1 + x3 + x5 + x6,
                        weight_method="ebal")
}
```
