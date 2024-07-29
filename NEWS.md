# adjustedCurves 0.9.0

* This is the first release of this package

# adjustedCurves 0.9.1

* Include CRAN in installation instructions
* Updated code and tests to run with updated versions of the mice and ggplot2 packages
* Updated documentation of package man page to include features not supported at the moment

# adjustedCurves 0.10.0

* Removed support for tmle, ostmle methods
* Changed citation information because manuscript was published
* Changed `print` method to be equal to `summary` method
* Fixed issues with unit-tests that require packages under "Suggests" only

# adjustedCurves 0.10.1

* Fixed small issues in unit tests caused by changes in the `WeightIt` package
* Made some small documentation updates

# adjustedCurves 0.11.0

Enhancements

* Added arguments `iso_reg` and `force_bounds` to `adjustedsurv()` and `adjustedcif()` functions to allow applying correction techniques outside plotting as well
* Added better support for multiple imputation when `variable`, `ev_time` or `event` contain missings, includes the new `mi_extrapolation` argument in `adjustedsurv()` and `adjustedcif()`
* Added the `ratio` argument to `adjusted_rmst()` and `adjusted_rmtl()` functions
* Added support for multiple `to` values in `adjusted_rmst()` and `adjusted_rmtl()`, which speeds up `plot_rmst_curve()` and `plot_rmtl_curve()` considerably

Bug Fixes

* There was a small bug in internal functions used to calculate integrals, which lead to slightly incorrect results whenever a survival curve reached 0 (or a CIF reached 1) and the "to" value was greater than the last observed time point. This may have impacted standard error estimates in `adjusted_rmst()`, `adjusted_rmtl()` and output of `adjusted_curve_test()` and functions that rely on those functions. This is fixed now, but may lead to slightly different estimates than in previous versions.

New Features

* Added the new methods `surv_tmle` and `cif_tmle`, based on the `concrete` package
* Added new instrumental variable based method `surv_iv_2SRIF`
* Added new methods `surv_prox_iptw`, `surv_prox_aiptw` based on code from Andrew Ying
* Added the `adjusted_curve_ratio()` function
* Added the `plot_curve_ratio()` function

Documentation

* Changed examples for the usage of `WeightIt` as suggested by Noah Greifer
* Added a new vignette with an overview of implemented features of each method
* Small changes to formulations

# adjustedCurves 0.11.1

Enhancements

* Re-factored internal code to vastly increase speed of bootstrapping related computations

New features

* Added risk table functionality for `plot.adjustedsurv()` (all arguments starting with `risk_table`)
* Allow estimation of difference and ratios in `plot_rmst_curve()` function
* Allow estimation of difference and ratios in `plot_rmtl_curve()` function
* Allow estimation of difference and ratios in `adjusted_surv_quantile()` function

Refactored

* Re-factored examples to only be executed if suggested packages are installed
* Renamed `adjsurv` and `adjcif` output objects of `adjustedsurv()` and `adjsutedcif()` respectively to `adj`
* Put functionality of `difference` and `ratio` arguments into one `contrast` argument in `adjusted_rmst()`, `adjusted_rmtl()`, `adjusted_surv_quantile()`, `plot_rmst_curve()` and `plot_rmtl_curve()` functions
* Temporarily removed support for `tmle` in `adjustedsurv()` and `adjustedcif()` due to `concrete` being removed from CRAN

Documentation

* Re-worked introduction vignette
* Added FAQ vignette
* Added Group Comparison vignette

# adjustedCurves 0.11.2

Bug Fixes

* Fixed a bug that resulted in additional arguments passed by the user through the three-dot (`...`) syntax not being correctly evaluated. For example, supplying `estimand="ATT"` to `adjustedsurv()` with `method="iptw_km"` and a formula in the `treatment_model` argument would result in the `estimand` argument not being passed to `weightit()`.
* Fixed issues that ocurred with `plot.adjustedsurv()` when the `adjustedsurv` object was created in a function or loop
* In previous versions the formula to pool standard errors when using multiple imputation was not implemented correctly. Fixed now, might lead to slightly different results when using multiply imputed data in the new and earlier versions. Many thanks to Dr. Jack M Wolf for finding and fixing this issue.
