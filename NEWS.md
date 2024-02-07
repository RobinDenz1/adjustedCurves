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

# adjustedCurves 0.10.2

* Changed examples for the usage of `WeightIt` as suggested by Noah Greifer
* Added the new methods `surv_tmle` and `cif_tmle`, based on the `concrete` package
* Added new instrumental variable based method `surv_iv_2SRIF`
* Added arguments `iso_reg` and `force_bounds` to `adjustedsurv` and `adjustedcif` functions
* Added new methods `surv_prox_iptw`, `surv_prox_aiptw` based on code from Andrew Ying
* Added a new vignette with an overview of implemented features of each method
* Added better support for multiple imputation when `variable`, `ev_time` or `event` contain missings
* Added the `ratio` argument to `adjusted_rmst()` and `adjusted_rmtl()` functions
* Added the `adjusted_curve_ratio()` function
* Added the `plot_curve_ratio()` function
