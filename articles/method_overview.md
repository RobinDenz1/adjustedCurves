# Overview of methods in adjustedCurves

The amount of methods implemented in this package can be overwhelming at
first, making one wonder which method should be used. This small
vignette exists to make this choice a little easier by providing a
concise overview of the functionality of each method implemented in the
[`adjustedsurv()`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md)
and
[`adjustedcif()`](https://robindenz1.github.io/adjustedCurves/reference/adjustedcif.md)
functions. Note that this vignette does not contain a description of how
these methods work or when. Information about that can be found in Denz
et al. (2023) or the respective documentation pages and the cited
literature therein.

## Methods in `adjustedsurv()`

The following table gives a general overview of all supported methods in
[`adjustedsurv()`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md):

|     | Method          | Supports Unmeasured Confounding | Supports Categorical Treatments | Supports Continuous Confounders | Approximate SE available | Always in Bounds | Always non-increasing | Doubly-Robust | Supports Dependent Censoring | Type of Adjustment | Is Nonparametric | Computation Speed | Dependencies     |
|:----|:----------------|:--------------------------------|:--------------------------------|:--------------------------------|:-------------------------|:-----------------|:----------------------|:--------------|:-----------------------------|:-------------------|:-----------------|:------------------|:-----------------|
| 1   | “direct”        | no                              | yes                             | yes                             | yes                      | yes              | yes                   | no            | no                           | outcome            | no               | \+                | riskRegression   |
| 2   | “direct_pseudo” | no                              | yes                             | yes                             | no                       | yes              | no                    | no            | yes                          | outcome            | no               | \- -              | geepack, prodlim |
| 3   | “iptw_km”       | no                              | yes                             | yes                             | yes                      | yes              | yes                   | no            | (no)                         | treatment          | depends          | ++                | \-               |
| 4   | “iptw_cox”      | no                              | yes                             | yes                             | no                       | yes              | yes                   | no            | (no)                         | treatment          | depends          | ++                | \-               |
| 5   | “iptw_pseudo”   | no                              | yes                             | yes                             | yes                      | no               | no                    | no            | yes                          | treatment          | depends          | \-                | prodlim          |
| 6   | “matching”      | no                              | no                              | yes                             | no                       | yes              | yes                   | no            | no                           | treatment          | depends          | \-                | Matching         |
| 7   | “emp_lik”       | no                              | no                              | yes                             | no                       | yes              | yes                   | no            | no                           | treatment          | yes              | \+                | MASS             |
| 8   | “aiptw”         | no                              | no                              | yes                             | yes                      | no               | no                    | yes           | yes                          | both               | no               | \-                | riskRegression   |
| 9   | “aiptw_pseudo”  | no                              | yes                             | yes                             | yes                      | no               | no                    | yes           | yes                          | both               | no               | \- -              | geepack, prodlim |
| 11  | “strat_amato”   | no                              | yes                             | no                              | no                       | yes              | yes                   | no            | no                           | \-                 | yes              | +++               | \-               |
| 12  | “strat_nieto”   | no                              | yes                             | no                              | yes                      | yes              | yes                   | no            | no                           | \-                 | yes              | +++               | \-               |
| 13  | “strat_cupples” | no                              | yes                             | no                              | no                       | yes              | yes                   | no            | no                           | \-                 | yes              | +++               | \-               |
| 14  | “iv_2SRIF”      | yes                             | no                              | yes                             | no                       | yes              | yes                   | no            | no                           | \-                 | no               | \+                | \-               |
| 15  | “prox_iptw”     | yes                             | no                              | yes                             | yes                      | no               | no                    | no            | no                           | treatment          | no               | \- -              | numDeriv         |
| 16  | “prox_aiptw”    | yes                             | no                              | yes                             | yes                      | no               | no                    | yes           | no                           | both               | no               | \- -              | numDeriv         |
| 17  | “km”            | no                              | yes                             | no                              | yes                      | yes              | yes                   | no            | no                           | none               | yes              | +++               | \-               |

For methods `"iptw_km"` and `"iptw_cox"` we wrote “(no)” in whether they
support dependent censoring, because there is no direct implementation
to handle it in this package. By supplying inverse probability of
censoring weights to the `treatment_model` argument it is, however,
possible to use those estimators to adjust for dependent censoring as
well. If both inverse probability of treatment (or more general
covariate balancing weights) **and** inverse probability of censoring
weights should be used, the user can simply multiply the subject-level
weights and supply the results to the `treatment_model` argument.

The following table gives an overview of the supported input to the
`treatment_model` argument for methods that require it:

| Method         | Allowed Input to treatment_model argument               |
|:---------------|:--------------------------------------------------------|
| “iptw_km”      | glm or multinom object, weights, formula for weightit() |
| “iptw_cox”     | glm or multinom object, weights, formula for weightit() |
| “iptw_pseudo”  | glm or multinom object, weights, formula for weightit() |
| “matching”     | glm object or propensity scores                         |
| “aiptw”        | glm object                                              |
| “aiptw_pseudo” | glm or multinom object or propensity scores             |

After having created an `adjustedsurv` object using the
[`adjustedsurv()`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md)
function, the following functions can be used to create plots, transform
the output or calculate further statistics:

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html): Plots the
  estimated adjusted survival curves
- [`adjusted_curve_diff()`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_curve_diff.md):
  Calculates differences in survival probabilities
- [`adjusted_curve_ratio()`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_curve_diff.md):
  Calculates ratios of survival probabilities
- [`plot_curve_diff()`](https://robindenz1.github.io/adjustedCurves/reference/plot_curve_diff.md):
  Plots differences in survival probabilities
- [`plot_curve_ratio()`](https://robindenz1.github.io/adjustedCurves/reference/plot_curve_diff.md):
  Plots ratios of survival probabilities
- [`adjusted_surv_quantile()`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_surv_quantile.md):
  Calculates adjusted survival time quantiles
- [`adjusted_rmst()`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_rmst.md):
  Calculates adjusted restricted mean survival times
- [`plot_rmst_curve()`](https://robindenz1.github.io/adjustedCurves/reference/plot_rmst_curve.md):
  Plots adjusted restricted mean survival time curves
- [`adjusted_rmtl()`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_rmtl.md):
  Calculates adjusted restricted mean time lost
- [`plot_rmtl_curve()`](https://robindenz1.github.io/adjustedCurves/reference/plot_rmtl_curve.md):
  Plots adjusted restricted mean time lost curves
- [`adjusted_curve_test()`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_curve_test.md):
  Performs a test of adjusted survival curve equality in an interval
- [`as_ggsurvplot_df()`](https://robindenz1.github.io/adjustedCurves/reference/as_ggsurvplot_df.md):
  Transforms the output to a concise `data.frame`

## Methods in `adjustedcif()`

The following table gives a general overview of all supported methods in
[`adjustedcif()`](https://robindenz1.github.io/adjustedCurves/reference/adjustedcif.md):

|     | Method           | Supports Unmeasured Confounding | Supports Categorical Treatments | Supports Continuous Confounders | Approximate SE available | Always in Bounds | Always non-increasing | Doubly-Robust | Supports Dependent Censoring | Type of Adjustment | Is Nonparametric | Computation Speed | Dependencies     |
|:----|:-----------------|:--------------------------------|:--------------------------------|:--------------------------------|:-------------------------|:-----------------|:----------------------|:--------------|:-----------------------------|:-------------------|:-----------------|:------------------|:-----------------|
| 1   | “direct”         | no                              | yes                             | yes                             | yes                      | yes              | yes                   | no            | no                           | outcome            | no               | \+                | riskRegression   |
| 2   | “direct_pseudo”  | no                              | yes                             | yes                             | no                       | yes              | no                    | no            | no                           | outcome            | no               | \- -              | geepack, prodlim |
| 3   | “iptw”           | no                              | yes                             | yes                             | yes                      | yes              | yes                   | no            | yes                          | treatment          | no               | \+                | riskRegression   |
| 4   | “iptw_pseudo”    | no                              | yes                             | yes                             | yes                      | no               | no                    | no            | no                           | treatment          | depends          | \+                | prodlim          |
| 5   | “matching”       | no                              | no                              | yes                             | no                       | yes              | yes                   | no            | no                           | treatment          | depends          | \-                | Matching         |
| 6   | “aiptw”          | no                              | no                              | yes                             | yes                      | no               | no                    | yes           | yes                          | both               | no               | \-                | riskRegression   |
| 7   | “aiptw_pseudo”   | no                              | yes                             | yes                             | yes                      | no               | no                    | yes           | no                           | both               | no               | \- -              | geepack, prodlim |
| 9   | “aalen_johansen” | no                              | yes                             | no                              | yes                      | yes              | yes                   | no            | no                           | none               | yes              | ++                | cmprsk           |

The following table gives an overview of the supported input to the
`treatment_model` argument for methods that require it:

| Method         | Allowed Input to treatment_model argument               |
|:---------------|:--------------------------------------------------------|
| “iptw”         | glm or multinom object                                  |
| “iptw_pseudo”  | glm or multinom object, weights, formula for weightit() |
| “matching”     | glm object or propensity scores                         |
| “aiptw”        | glm object                                              |
| “aiptw_pseudo” | glm or multinom object or propensity scores             |

Note that method `"iptw"` currently does not support directly supplying
weights or propensity scores. This is due to it relying on the `ate`
function of the `riskRegression` package, which only accepts glm or
multinom objects. This may be changed in the future.

After having created an `adjustedcif` object using the
[`adjustedcif()`](https://robindenz1.github.io/adjustedCurves/reference/adjustedcif.md)
function, the following functions can be used to create plots, transform
the output or calculate further statistics:

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html): Plots the
  estimated adjusted CIFs
- [`adjusted_curve_diff()`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_curve_diff.md):
  Calculates differences in CIFs
- [`adjusted_curve_ratio()`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_curve_diff.md):
  Calculates ratios of CIFs
- [`plot_curve_diff()`](https://robindenz1.github.io/adjustedCurves/reference/plot_curve_diff.md):
  Plots differences in CIFs over time
- [`plot_curve_ratio()`](https://robindenz1.github.io/adjustedCurves/reference/plot_curve_diff.md):
  Plots ratios of survival probabilities
- [`adjusted_rmtl()`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_rmtl.md):
  Calculates adjusted restricted mean time lost
- [`plot_rmtl_curve()`](https://robindenz1.github.io/adjustedCurves/reference/plot_rmtl_curve.md):
  Plots adjusted restricted mean time lost curves
- [`adjusted_curve_test()`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_curve_test.md):
  Performs a test of adjusted CIF equality in an interval

## References

Robin Denz, Renate Klaaßen-Mielke, and Nina Timmesfeld (2023). “A
Comparison of Different Methods to Adjust Survival Curves for
Confounders”. In: Statistics in Medicine 42.10, pp. 1461-1479
