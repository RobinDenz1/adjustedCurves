# Augmented Inverse Probability of Treatment Weighted CIFs

This page explains the details of estimating augmented inverse
probability of treatment weighted cumulative incidence functions for
competing risks data (`method="aiptw"` in the
[`adjustedcif`](https://robindenz1.github.io/adjustedCurves/reference/adjustedcif.md)
function). All regular arguments of the `adjustedcif` function can be
used. Additionally, the `outcome_model` argument and the
`treatment_model` argument have to be specified in the `adjustedcif`
call. Further arguments specific to this method are listed below.

## Arguments

- outcome_model:

  \[**required**\] Must be a `CauseSpecificCox` model object created
  using the [`CSC`](https://rdrr.io/pkg/riskRegression/man/CSC.html)
  function, modeling the time-to-event mechanism. See details and
  examples.

- treatment_model:

  \[**required**\] Must be a `glm` model object with `variable` as
  response variable. See details and examples.

- censoring_model:

  Must be a `coxph` model object, modeling the censoring mechanism or
  `NULL`. If `NULL` (default) independent censoring is assumed. See
  details and examples.

- verbose:

  Whether to print estimation information of the `ate` function in the
  riskRegression package. Defaults to `FALSE`.

- ...:

  Further arguments passed to
  [`ate`](https://rdrr.io/pkg/riskRegression/man/ate.html).

## Details

- **Type of Adjustment:** Requires both a treatment assignment model
  ([`glm`](https://rdrr.io/r/stats/glm.html)) and a outcome model
  ([`CSC`](https://rdrr.io/pkg/riskRegression/man/CSC.html)). Also
  allows, but does not rely on, an additional model describing the
  censoring mechanism (a
  [`coxph`](https://rdrr.io/pkg/survival/man/coxph.html) object).

- **Doubly-Robust:** Estimates are Doubly-Robust.

- **Categorical groups:** Currently only two groups in `variable` are
  allowed. Must still be a factor variable.

- **Approximate Variance:** Calculations to approximate the variance and
  confidence intervals are available.

- **Allowed Time Values:** Allows both continuous and integer time.

- **Bounded Estimates:** Estimates are not guaranteed to be bounded in
  the 0 to 1 probability range.

- **Monotone Function:** Estimates are not guaranteed to be monotone.

- **Dependencies:** This method relies on the riskRegression package.

Instead of only modeling the outcome mechanism or the treatment
assignment mechanism, both kind of models are required to use this
method. If either of those models are correctly specified, unbiased
estimates will be obtained. Can also be used to adjust for dependent
censoring using a Cox-Regression model. An obvious advantage of this
method is it's doubly robust property. This however comes at the price
of some efficiency. It is also possible that some estimates fall outside
the 0 and 1 probability bounds, particularly if the time is near 0 or
the maximal observed event time. There is also no guarantee that the
estimated CIFs will be monotonically increasing. For more information on
the methods the user is referred to the literature listed in the
references.

This function is basically just a wrapper around the
[`ate`](https://rdrr.io/pkg/riskRegression/man/ate.html) function from
the riskRegression package. Additional arguments may be passed to that
function using the `...` syntax. It is however recommended to use
[`ate`](https://rdrr.io/pkg/riskRegression/man/ate.html) directly in
these cases.

## Value

Adds the following additional objects to the output of the `adjustedcif`
function:

- `ate_object`: The object returned by the `ate` function.

## References

James M. Robins and Andrea Rotnitzky (1992). "Recovery of Information
and Adjustment for Dependent Censoring Using Surrogate Markers". In:
AIDS Epidemiology: Methodological Issues. Ed. by Nicholas P. Jewell,
Klaus Dietz, and Vernon T. Farewell. New York: Springer Science +
Business Media, pp. 297-331

Alan E. Hubbard, Mark J. van der Laan, and James M. Robins (2000).
"Nonparametric Locally Efficient Estimation of the Treatment Specific
Survival Distribution with Right Censored Data and Covariates in
Observational Studies". In: Statistical Models in Epidemiology, the
Environment, and Clinical Trials. Ed. by M. Elizabeth Halloran and
Donald Berry. New York: Springer Science + Business Media, pp. 135-177

Brice Maxime Hugues Ozenne, Thomas Harder Scheike, and Laila Staerk
(2020). "On the Estimation of Average Treatment Effects with
Right-Censored Time to Event Outcome and Competing Risks". In:
Biometrical Journal 62, pp. 751-763

## Author

The wrapper function was written by Robin Denz, the `ate` function
(which this wrapper is build around) was written by other people. See
[`?ate`](https://rdrr.io/pkg/riskRegression/man/ate.html) for more
details.

## See also

[`ate`](https://rdrr.io/pkg/riskRegression/man/ate.html),
[`CSC`](https://rdrr.io/pkg/riskRegression/man/CSC.html),
[`coxph`](https://rdrr.io/pkg/survival/man/coxph.html),
[`glm`](https://rdrr.io/r/stats/glm.html)

## Examples

``` r
library(adjustedCurves)
library(survival)

if (requireNamespace("riskRegression")) {

library(riskRegression)

set.seed(42)

# simulate some data as example
sim_dat <- sim_confounded_crisk(n=50, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# estimate a cause-specific cox-regression for the outcome
cox_mod <- CSC(Hist(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6 + group,
               data=sim_dat)

# estimate a treatment assignment model
glm_mod <- glm(group ~ x1 + x3 + x5 + x6, data=sim_dat, family="binomial")

# use it to calculate adjusted survival curves
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      cause=1,
                      method="aiptw",
                      outcome_model=cox_mod,
                      treatment_model=glm_mod,
                      conf_int=FALSE)

# plot the curves
plot(adjcif)
}
#> Warning: Rare event 
#> Ignoring unknown labels:
#> • linetype : "Group"
#> • fill : "Group"
```
