# Using Propensity-Score Matching to Calculate Adjusted Survival Curves

This page explains the details of estimating adjusted survival curves
using propensity-score matching for single event time-to-event data
(`method="matching"` in the
[`adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md)
function). All regular arguments of the `adjustedsurv` function can be
used. Additionally, the `treatment_model` argument has to be specified
in the `adjustedsurv` call. Further arguments specific to this method
are listed below.

## Arguments

- treatment_model:

  \[**required**\] Must be either a model object with `variable` as
  response variable or a vector of previously estimated propensity
  scores.

- gtol:

  Tolerance at which estimated treatment assignment probabilities are
  truncated. Every propensity score bigger than 1 - `gtol` is set to 1 -
  `gtol` and every propensity score smaller than `gtol` is set to
  `gtol`. Useful when there are extreme propensity scores close to 0
  or 1. Defaults to 0.001.

- ...:

  Further arguments passed to the
  [`Match`](https://rdrr.io/pkg/Matching/man/Match.html) function of the
  Matching Package.

## Details

- **Type of Adjustment:** Requires a model describing the treatment
  assignment mechanism. This must be either a
  [`glm`](https://rdrr.io/r/stats/glm.html) object.

- **Doubly-Robust:** Estimates are not Doubly-Robust.

- **Categorical groups:** Only two groups in `variable` are allowed.
  Must be a factor variable with exactly two levels.

- **Approximate Variance:** Calculations to approximate the variance and
  confidence intervals are currently not available. Bootstrap confidence
  intervals can however be calculated with all supported models. See
  [`?adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md)
  for more information on bootstrapping.

- **Allowed Time Values:** Allows both continuous and integer time.

- **Bounded Estimates:** Estimates are guaranteed to be bounded in the 0
  to 1 probability range.

- **Monotone Function:** Estimates are guaranteed to be monotone.

- **Dependencies:** This method relies on the Matching package.

Using the estimated propensity score, the individual observations in the
dataset are matched to each other creating a new dataset in which the
covariate distributions are balanced in respect to the two groups
defined by `variable`. A simple Kaplan-Meier estimator is then used to
calculate the confounder-adjusted survival curves. This corresponds to
the method described in Austin (2014). Details on the algorithm used for
matching can be found in the documentation of the Matching package.

We choose not to implement other matching based estimators (see Winnett
& Sasieni (2002), Galimberti et al. (2002) and Austin (2020)) because of
the wide range of matching algorithms and parameters. Trying to automate
the matching process in a function like this would, in our opinion,
disrupt the workflow of the user while also encouraging suboptimal
practices. We however included this simple version of a matching
estimator as a reference and to raise the awareness that using matching
is a valid method to obtain adjusted survival curves.

Simulation studies have shown that this particular method as implemented
here is significantly less efficient than other methods included in this
R-Package. While it does produce unbiased estimates, the variation in
these estimates is very high. We suggest using one of the other
available methods.

## Value

Adds the following additional objects to the output of the
`adjustedsurv` function:

- `match_object`: The object creates using the `Match` function.

- `survfit_object`: The `survfit` object fit on the matched data.

## References

Peter C. Austin (2014). "The Use of Propensity Score Methods with
Survival or Time-To-Event Outcomes: Reporting Measures of Effect Similar
to those Used in Randomized Experiments". In: Statistics in Medicine 33,
pp. 1242-1258

Angela Winnett and Peter Sasieni (2002). "Adjusted Nelson-Aalen
Estimates with Retrospective Matching". In: Journal of the American
Statistical Association 97.457, pp. 245-256

Stefania Galimberti, Peter Sasieni, and Maria Grazia Valsecchi (2002).
"A Weighted Kaplan-Meier Estimator for Matched Data with Application to
the Comparison of Chemotherapy and Bone-Marrow Transplant in Leukaemia".
In: Statistics in Medicine 21, pp. 3847-3864

Peter C. Austin, Neal Thomas, and Donald B. Rubin (2020).
"Covariate-Adjusted Survival Analyses in Propensity-Score Matched
Samples: Imputing Potential Time- To-Event Outcomes". In: Statistical
Methods in Medical Research 29.3, pp. 728-751

## Author

Robin Denz

## See also

[`Match`](https://rdrr.io/pkg/Matching/man/Match.html),
[`survfit`](https://rdrr.io/pkg/survival/man/survfit.html)

## Examples

``` r
library(adjustedCurves)

if (requireNamespace("Matching")) {

library(Matching)

set.seed(42)

# simulate some data as example
sim_dat <- sim_confounded_surv(n=50, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# estimate treatment assignment model
glm_mod <- glm(group ~ x1 + x2 + x4 + x6, data=sim_dat, family="binomial")

# calculate adjusted survival curves
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="matching",
                        treatment_model=glm_mod)

# Alternatively, supply the propensity score directly
# Here we use the logistic regression to calculate it, so we get
# exactly the same result. The propensity score can be calculated in
# any other way in practice, allowing flexibility
ps_score <- glm_mod$fitted.values

adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="matching",
                        treatment_model=ps_score)

# plot the curves
plot(adjsurv)
}
#> Ignoring unknown labels:
#> • linetype : "Group"
#> • fill : "Group"
```
