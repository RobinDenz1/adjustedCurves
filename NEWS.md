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

* Added the `difference` and `ratio` arguments to `plot_rmst_curve()` and `plot_rmtl_curve()`
* Added the `difference` and `ratio` functionality to the `adjusted_surv_quantile()` function
* Re-factored internal code to vastly increase speed of bootstrapping related computations
* Added risk table functionality for `plot.adjustedsurv()`
* Renamed `adjsurv` and `adjcif` output objects of `adjustedsurv()` and `adjsutedcif()` respectively to `adj`
