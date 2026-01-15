# Estimate Cause-Specific Confounder-Adjusted Cumulative Incidence Functions

This is one of two main functions of this R-Package. It allows the user
to estimate cause-specific confounder-adjusted cumulative incidence
functions in the presence of competing events using a variety of
different methods. Some of these methods require additional packages to
be installed and, depending on the specified method, there might be
additional required arguments in the function call. More information is
available on the documentation page of the respective `cif_method`.

## Usage

``` r
adjustedcif(data, variable, ev_time, event, cause,
            method, conf_int=FALSE, conf_level=0.95,
            times=NULL, bootstrap=FALSE, n_boot=500,
            n_cores=1, na.action=options()$na.action,
            clean_data=TRUE, iso_reg=FALSE,
            force_bounds=FALSE, mi_extrapolation=FALSE,
            ...)
```

## Arguments

- data:

  A `data.frame` object containing the needed time-to-event data in
  standard format. Can also be a `mids` object created with the mice
  package. See details for how this works.

- variable:

  A character string specifying the variable by which the cumulative
  incidence functions should be grouped. Must be a valid column name of
  `data`. The variable specified should needs to be a factor variable.

- ev_time:

  A character string specifying the variable indicating the
  time-to-event or time-to-censoring. Must be a valid column name of
  `data`.

- event:

  A character string specifying the numeric event indicator. The
  censoring indicator should be coded as 0 and all other events of
  interest as 1, 2, etc. Must be a valid column name of `data`.

- cause:

  The cause of interest for which the cumulative incidence functions
  should be estimated. Should be a number that appears in the `event`
  column of `data`.

- method:

  A character string specifying the adjustment method to use. Case
  sensitive. See details.

- conf_int:

  A logical variable, indicating whether the asymptotic variances and
  confidence intervals of the cumulative incidence should be estimated.
  Not available for all methods. More information can be found in the
  documentation of each method. For an alternative way to get confidence
  intervals, see the `bootstrap` argument.

- conf_level:

  A number specifying the confidence level of asymptotic and/or
  bootstrap confidence intervals.

- times:

  A numeric vector of time points at which the cumulative incidences
  should be estimated or `NULL`. If `NULL` the cumulative incidence is
  estimated at all points in time at which any event occurred in the
  pooled sample.

- bootstrap:

  A logical variable indicating whether bootstrapping should be
  performed or not. In bootstrapping, a number of simple random samples
  with replacement of size `nrow(data)` are drawn from `data`. For each
  sample the calculations are repeated and used to estimate standard
  errors and confidence intervals. This can be used to obtain confidence
  intervals when asymptotic variance calculations are not available.

- n_boot:

  Number of bootstrap replications to perform. Ignored if `bootstrap` is
  `FALSE`.

- n_cores:

  The number of cores to use when calculating bootstrap estimates.
  Ignored if `bootstrap=FALSE`. Is set to 1 by default, resulting in
  single threaded processing. Internally uses the doParallel package if
  `n_cores > 1`. In that case it also uses the doRNG package to make the
  results replicable. See
  [`?doRNG`](https://rdrr.io/pkg/doRNG/man/infoDoRNG.html) and
  `?doParallel` for more details. Using multiple cores will speed up the
  calculation considerably in most cases.

- na.action:

  How missing values should be handled. Can be one of: na.fail, na.omit,
  na.pass or na.exclude. Also accepts strings of the function names. See
  [`?na.action`](https://rdrr.io/r/stats/na.action.html) for more
  details. By default it uses the na.action which is set in the global
  options by the respective user.

- clean_data:

  If `TRUE` all columns which are not needed for the estimation are
  removed from `data` before any further calculations are performed.
  This ensures that calls to `na.omit` (see argument `na.action`) do not
  remove rows which are fully observed in respect to relevant columns
  due to missing values in irrelevant columns. Set to `FALSE` to skip
  this step. Usually this argument can be ignored. When using
  non-standard outcome models however it should be set to `FALSE`.

- iso_reg:

  Either `TRUE` or `FALSE` (default), controlling whether isotonic
  regression is performed on the resulting failure probability
  estimates. This can be used to ensure that the CIFs are
  non-decreasing. Since only a few methods may exhibit this problem,
  this argument is only relevant for some methods (see method specific
  documentation).

- force_bounds:

  Either `TRUE` or `FALSE` (default), controlling whether the resulting
  failure probability estimates should be forced to lie between 0 and 1.
  If `TRUE` and there are values higher than 1, they are simply set
  to 1. Values lower than 0 are similarly set to 0. Since only a few
  methods may exhibit this problem, this argument is only relevant for
  some methods (see method specific documentation).

- mi_extrapolation:

  Whether to allow extrapolation due to imputed survival times or not.
  This argument is only relevant when using multiply imputed `data` with
  missing covariates in `variable`, `ev_time` or `event`. Depending on
  the algorithm used to obtain the imputed datasets, it may be possible
  that one or more imputed datasets contain survival times in a group
  that are larger than the maximum observed survival time in that group.
  This may lead to unwanted extrapolation (e.g. the survival curves
  extending further than they should). By keeping this argument at
  `FALSE`, these times are removed from the output. If set to `TRUE`,
  all available estimates will be used.

- ...:

  Further arguments passed to the respective `cif_method`. For example
  when using `method="direct"` all further arguments are passed to the
  `cif_direct` function. See details.

## Details

The primary purpose of the `adjustedcif` function is to provide a
convenient way to estimate confounder-adjusted cumulative incidence
functions using any of the methods provided in the literature. A
[`plot`](https://robindenz1.github.io/adjustedCurves/reference/plot.adjustedcif.md)
method is provided to graphically display the estimated cumulative
incidence functions as well. Currently the following methods can be
used:

- "[direct](https://robindenz1.github.io/adjustedCurves/reference/cif_direct.md)":
  Direct Standardization based on a model describing the outcome
  mechanism ([`CSC`](https://rdrr.io/pkg/riskRegression/man/CSC.html),
  [`FGR`](https://rdrr.io/pkg/riskRegression/man/FGR.html), ..).

- "[direct_pseudo](https://robindenz1.github.io/adjustedCurves/reference/cif_direct_pseudo.md)":
  Direct Standardization based on Pseudo-Values.

- "[iptw](https://robindenz1.github.io/adjustedCurves/reference/cif_iptw.md)":
  A weighted Aalen-Johansen estimator.

- "[iptw_pseudo](https://robindenz1.github.io/adjustedCurves/reference/cif_iptw_pseudo.md)":
  A weighted estimator based on Pseudo-Values.

- "[matching](https://robindenz1.github.io/adjustedCurves/reference/cif_matching.md)":
  Using Propensity Score Matching to estimate the adjusted CIFs.

- "[aiptw](https://robindenz1.github.io/adjustedCurves/reference/cif_aiptw.md)":
  An Augmented Inverse Probability of Treatment Weighting estimator.

- "[aiptw_pseudo](https://robindenz1.github.io/adjustedCurves/reference/cif_aiptw_pseudo.md)":
  An Augmented Inverse Probability of Treatment Weighting estimator
  using Pseudo-Values.

- "[aalen_johansen](https://robindenz1.github.io/adjustedCurves/reference/cif_aalen_johansen.md)":
  A simple stratified Aalen-Johansen estimator without any form of
  adjustment.

A short description of each method is contained in the documentation of
the respective `cif_method` function. A concise overview of the
supported functionality of each method can be found in the associated
vignette
([`vignette(topic="method_overview", package="adjustedCurves")`](https://robindenz1.github.io/adjustedCurves/articles/method_overview.md)).
For more detailed descriptions the cited literature in the respective
documentation pages can be used. The documentation for `method="direct"`
for example can be accessed using
[`?cif_direct`](https://robindenz1.github.io/adjustedCurves/reference/cif_direct.md).

***Required & Optional Arguments***

Every method requires the specification of the `data`, `variable`,
`ev_time`, `event`, `cause` and `method` arguments. All other arguments
mentioned on this page are optional and work for all methods. Depending
on the method used, other arguments are required as well. Those can be
found on the top of the help page of the respective method. The help
pages also list additional optional arguments.

***Confidence Intervals***

For most methods approximations for the asymptotic variance of point
estimates of the CIF have been proposed in the literature. Where
available, those can be estimated and added to the output object using
`conf_int=TRUE`. It is however recommended to use bootstrapping to
estimate the variance instead, which can be done by setting
`bootstrap=TRUE`. The `n_boot` argument is set to 500 by default. This
number was chosen because it worked well in simulations but it does not
guarantee convergence in practice. Users are recommended to inspect the
bootstrapped estimates and adjust the number of replications
accordingly. To allow faster bootstrapping the user can choose to run
the function on multiple CPU cores in parallel using the `n_cores`
argument.

***Missing Data***

There are two ways to deal with missing data using this function. The
first is using the `na.action` argument. It simply calls the respective
`na.action` function on the data before doing any further processing. By
using `na.action="na.omit"` for example, only rows with complete data
are kept for the analysis.

Alternatively, this function also supports the use of multiple
imputation via the mice package. Instead of supplying a single
data.frame, the user should create a `mids` object using the `mice`
function and directly pass this to the `data` argument. When methods are
used which rely on previously estimated treatment or outcome models such
as `"direct"` or `"aiptw"`, the user is required to supply a `mira`
object instead of a single model. In other words: the models have to be
fit on every imputed dataset before supplying them to this function. See
[`?mice`](https://amices.org/mice/reference/mice.html) and the
associated documentation for more information on how to use multiple
imputation. When using `bootstrap=TRUE` and multiple imputation, the
bootstrapping is performed on every imputed dataset separately.
Cumulative Incidences are simply averaged across the imputed datasets
according to Rubins Rule. Confidence intervals are calculated by first
averaging the standard errors over all imputed datasets and afterwards
using this pooled value to obtain a new confidence interval with the
normal approximation.

***Competing Risks***

This function is meant to be used for data containing multiple competing
risks. If the data does not contain competing-events, it is recommended
to use the `adjustedsurv` function instead. It does not estimate the CIF
directly, but the CIF can be calculated from the survival using CIF =
1 - \\S(t)\\. This can be done automatically in the `plot.adjustedsurv`
function using `cif=TRUE`.

***Graphical Displays***

A general plot of the estimated adjusted CIFs can be obtained using the
associated
[`plot`](https://robindenz1.github.io/adjustedCurves/reference/plot.adjustedcif.md)
method. In addition, a plot of the difference between two estimated
adjusted CIFs can be produced using the
[`plot_curve_diff`](https://robindenz1.github.io/adjustedCurves/reference/plot_curve_diff.md)
function.

## Value

Returns an `adjustedcif` object containing the following objects:

- adj:

  A `data.frame` of estimated cumulative incidences for `cause` at some
  points in time for each level of `variable`. Depending on the
  arguments used also includes standard errors and confidence intervals.

- method:

  The method used to adjust the CIFs.

- categorical:

  Whether there are more than 2 groups in `variable`.

- call:

  The original function call.

When the argument `bootstrap` is set to `TRUE` it additionally contains
the following objects:

- boot_data:

  The adjusted CIFs estimated in each bootstrap sample.

- boot_adj:

  The mean CIFs of all bootstrap samples and corresponding standard
  errors and percentile confidence intervals.

When multiple imputation was used, the function additionally contains a
`mids_analyses` object, containing the `adjustedcif` objects for each
imputed dataset.

Some method specific objects might also be contained in the output.

## References

Robin Denz, Renate Klaaßen-Mielke, and Nina Timmesfeld (2023). "A
Comparison of Different Methods to Adjust Survival Curves for
Confounders". In: Statistics in Medicine 42.10, pp. 1461-1479

More relevant literature can be found in the respective `cif_method`
documentation.

## Author

The function itself was written by Robin Denz, but some `cif_method`
functions include wrappers for functions written by other people. More
information can be found in the respective `cif_method` documentation.

## See also

[`plot.adjustedcif`](https://robindenz1.github.io/adjustedCurves/reference/plot.adjustedcif.md),
[`adjusted_rmtl`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_rmtl.md),
[`plot_curve_diff`](https://robindenz1.github.io/adjustedCurves/reference/plot_curve_diff.md)

## Examples

``` r
library(adjustedCurves)
library(survival)

set.seed(42)

# simulate some example data
sim_dat <- sim_confounded_crisk(n=50, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# treatment assignment model
glm_mod <- glm(group ~ x2 + x3 + x5 + x6, data=sim_dat, family="binomial")

if (requireNamespace("riskRegression")) {

library(riskRegression)

# outcome model
cox_mod <- CSC(Hist(time, event) ~ x1 + x2 + x4 + x5 + group, data=sim_dat)

# using direct adjustment with asymptotic confidence intervals for cause 1
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      cause=1,
                      method="direct",
                      outcome_model=cox_mod,
                      conf_int=TRUE,
                      bootstrap=FALSE)

# using IPTW with asymptotic confidence intervals for cause 2
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="iptw",
                      cause=2,
                      treatment_model=glm_mod,
                      conf_int=TRUE,
                      bootstrap=FALSE)

# using AIPTW with asymptotic confidence intervals for cause 1
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      cause=1,
                      method="aiptw",
                      outcome_model=cox_mod,
                      treatment_model=glm_mod,
                      conf_int=TRUE,
                      bootstrap=FALSE)

# using direct adjustment at custom points in time
custom_times <- c(0.001, 0.1, 0.2, 0.6, 1.1)
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      cause=1,
                      method="direct",
                      outcome_model=cox_mod,
                      conf_int=TRUE,
                      bootstrap=FALSE,
                      times=custom_times)

# using bootstrapping with direct adjustment
# NOTE: In practice the number of bootstrap replications should be
#       greater than 10. This is only shown here for convenience.
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      cause=1,
                      method="direct",
                      outcome_model=cox_mod,
                      conf_int=TRUE,
                      bootstrap=TRUE,
                      n_boot=10)

}
#> Warning: Rare event 
#> Warning: Rare event 
#> Warning: Rare event 
#> Warning: Rare event 
#> Warning: Loglik converged before variable  5 ; coefficient may be infinite. 
#> Warning: Rare event 
#> Warning: Estimated risk outside the range [0,1].
#> Consider setting the argument 'product.limit' to FALSE. 
#> Warning: Loglik converged before variable  5 ; coefficient may be infinite. 
#> Warning: Rare event 
#> Warning: Rare event 
#> Warning: Estimated risk outside the range [0,1].
#> Consider setting the argument 'product.limit' to FALSE. 
#> Warning: Rare event 
#> Warning: Estimated risk outside the range [0,1].
#> Consider setting the argument 'product.limit' to FALSE. 
#> Warning: Loglik converged before variable  5 ; coefficient may be infinite. 
#> Warning: Rare event 
#> Warning: Rare event 
#> Warning: Rare event 
#> Warning: Estimated risk outside the range [0,1].
#> Consider setting the argument 'product.limit' to FALSE. 
#> Warning: Rare event 
#> Warning: Rare event 
#> Warning: Estimated risk outside the range [0,1].
#> Consider setting the argument 'product.limit' to FALSE. 
#> Warning: Rare event 

# not run because those are too slow
# \donttest{
# using bootstrapping with direct adjustment, run in parallel
# on two cores
library(foreach)
library(parallel)
library(doRNG)
#> Loading required package: rngtools

if (requireNamespace("riskRegression")) {

adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      cause=1,
                      method="direct",
                      outcome_model=cox_mod,
                      conf_int=TRUE,
                      bootstrap=TRUE,
                      n_boot=4,
                      n_cores=2)
}
#> Loading required namespace: doParallel
#> Warning: Rare event 

if (requireNamespace("mice") & requireNamespace("WeightIt")) {

# using multiple imputation
library(mice)
library(WeightIt)

# simulate some data as example
sim_dat <- sim_confounded_crisk(n=100, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# introduce random missingness in x1 as example
# NOTE: This is only done as an example, in reality you would
#       already have missing data, not introduce it yourself.
sim_dat$x1 <- ifelse(runif(n=50) < 0.5, sim_dat$x1, NA)

# perform multiple imputation
mids <- mice::mice(data=sim_dat, method="pmm", m=2, printFlag=FALSE)

# IPTW Pseudo using WeightIt on imputed data, for cause = 1
adj <- adjustedcif(data=mids,
                   variable="group",
                   ev_time="time",
                   event="event",
                   method="iptw_pseudo",
                   cause=1,
                   treatment_model=group ~ x1 + x2 + x5 + x6,
                   weight_method="ps")
plot(adj, force_bounds=TRUE, iso_reg=TRUE)
}
#> Loading required namespace: WeightIt
#> Ignoring unknown labels:
#> • linetype : "Group"
#> • fill : "Group"


# More specific examples can be found in the documentation of each
# respective cif_method. See ?cif_ + "method" for more information.
# }
```
