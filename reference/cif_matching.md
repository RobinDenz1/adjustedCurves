# Using Propensity-Score Matching to Calculate Adjusted CIFs

This page explains the details of estimating adjusted cumulative
incidence functions using propensity-score matching in a competing risks
setting (`method="matching"` in the
[`adjustedcif`](https://robindenz1.github.io/adjustedCurves/reference/adjustedcif.md)
function). All regular arguments of the `adjustedcif` function can be
used. Additionally, the `treatment_model` argument has to be specified
in the `adjustedcif` call. Further arguments specific to this method are
listed below.

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
  or 1. Defaults to 0.001,

- ...:

  Further arguments passed to the `Match` function of the Matching
  Package.

## Details

- **Type of Adjustment:** Requires a model describing the treatment
  assignment mechanism. This must be either a
  [`glm`](https://rdrr.io/r/stats/glm.html) object or a vector of
  propensity scores.

- **Doubly-Robust:** Estimates are not Doubly-Robust.

- **Categorical groups:** Only two groups in `variable` are allowed.
  Must be a factor variable with exactly two levels.

- **Approximate Variance:** Calculations to approximate the variance and
  confidence intervals are currently not available. Bootstrapping can
  still be used to estimate the confidence intervals (see
  [`?adjustedcif`](https://robindenz1.github.io/adjustedCurves/reference/adjustedcif.md)).

- **Allowed Time Values:** Allows both continuous and integer time.

- **Bounded Estimates:** Estimates are guaranteed to be bounded in the 0
  to 1 probability range.

- **Monotone Function:** Estimates are guaranteed to be monotone.

- **Dependencies:** This method relies on both the Matching and the
  cmprsk packages.

Using the estimated propensity score, the individual observations in the
dataset are matched to each other creating a new dataset in which the
covariate distributions are balanced in respect to the two groups
defined by `variable`. A simple Aalen-Johansen estimator is then used to
calculate the confounder-adjusted CIFs. This corresponds to the method
described in Austin & Fine (2019). Details on the algorithm used for
matching can be found in the documentation of the Matching package.

Simulation results showed that this specific implementation of this
method is the least efficient method contained in this R-Package. While
it does produce unbiased estimates, the variation in these estimates is
very high. We strongly suggest using one of the other methods
implemented here.

## Value

Adds the following additional objects to the output of the `adjustedcif`
function:

- `match_object`: The object creates using the `Match` function.

- `cuminc_object`: The `cuminc` object fit on the matched data.

## References

Peter C. Austin and Jason P. Fine (2019). "Propensity-Score Matching
with Competing Risks in Survival Analysis". In: Statistics in Medicine
38, pp. 751-777

## Author

Robin Denz

## See also

[`Match`](https://rdrr.io/pkg/Matching/man/Match.html),
[`cuminc`](https://rdrr.io/pkg/cmprsk/man/cuminc.html)

## Examples

``` r
library(adjustedCurves)
library(survival)

if (requireNamespace("Matching")) {

library(Matching)

set.seed(42)

# simulate some data as example
sim_dat <- sim_confounded_crisk(n=50, max_t=5)
sim_dat$group <- as.factor(sim_dat$group)

# estimate treatment assignment model
glm_mod <- glm(group ~ x1 + x2 + x4 + x6, data=sim_dat, family="binomial")

# calculate adjusted CIFs
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      cause=1,
                      method="matching",
                      treatment_model=glm_mod)
plot(adjcif)

# Alternatively, supply the propensity score directly
# Here we use the logistic regression to calculate it, so we get
# exactly the same result. The propensity score can be calculated in
# any other way in practice, allowing flexibility
ps_score <- glm_mod$fitted.values

adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      cause=1,
                      method="matching",
                      treatment_model=ps_score)

# plot the curves
plot(adjcif)
}
#> Loading required namespace: Matching
#> Loading required package: MASS
#> ## 
#> ##  Matching (Version 4.10-15, Build Date: 2024-10-14)
#> ##  See https://www.jsekhon.com for additional documentation.
#> ##  Please cite software as:
#> ##   Jasjeet S. Sekhon. 2011. ``Multivariate and Propensity Score Matching
#> ##   Software with Automated Balance Optimization: The Matching package for R.''
#> ##   Journal of Statistical Software, 42(7): 1-52. 
#> ##
#> Ignoring unknown labels:
#> • linetype : "Group"
#> • fill : "Group"
```
