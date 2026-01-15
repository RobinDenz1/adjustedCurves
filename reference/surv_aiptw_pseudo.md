# Augmented Inverse Probability of Treatment Weighted Survival Curves using Pseudo-Values

This page explains the details of estimating augmented inverse
probability of treatment weighted survival curves using Pseudo-Values
for single event time-to-event data (`method="aiptw_pseudo"` in the
[`adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md)
function). All regular arguments of the `adjustedsurv` function can be
used. Additionally, the `outcome_vars` argument and the
`treatment_model` argument have to be specified in the `adjustedsurv`
call. Further arguments specific to this method are listed below.

## Arguments

- outcome_vars:

  \[**required**\] A character vector of column names specifying
  variables to be used when modeling the outcome mechanism using
  [`geese`](https://rdrr.io/pkg/geepack/man/geese.html). See details and
  examples.

- treatment_model:

  \[**required**\] Must be a `glm` or `multinom` model object with
  `variable` as response variable. Alternatively you can supply a
  numeric vector of propensity scores directly. See details and
  examples.

- type_time:

  A character string specifying how the time should be modeled. Possible
  values are `"factor"` (modeling each point in time as a separate
  variable, the default), `"bs"` (modeling time using B-Splines) or
  `"ns"` (modeling time using natural splines).

- spline_df:

  The number of degrees of freedom used for the natural-spline or
  B-spline function. Ignored if `type_time="factor"`. Defaults to 5.

- censoring_vars:

  An optional character vector specifying variables in `data`. Those are
  used in the calculation of inverse probability of censoring weighted
  pseudo observations. See `?pseudo_aareg` for more information. Set to
  `NULL` (default) to use standard pseudo-values without corrections for
  dependent censoring instead.

- ipcw_method:

  The specific method used in the calculation of inverse probability of
  censoring weighted pseudo observations. Can be either `"binder"`
  (default) or `"hajek"`. See `?pseudo_aareg` for more information.
  Ignored if `censoring_vars=NULL`.

## Details

- **Type of Adjustment:** Requires a treatment assignment model
  ([`glm`](https://rdrr.io/r/stats/glm.html) or
  [`multinom`](https://rdrr.io/pkg/nnet/man/multinom.html)) and a
  character vector of variable names used to model the outcome mechanism
  (internally uses
  [`geese`](https://rdrr.io/pkg/geepack/man/geese.html)).
  Covariate-Dependent censoring can be corrected for using inverse
  probability of censoring weighted pseudo-values (Binder et al. 2014)

- **Doubly-Robust:** Estimates are Doubly-Robust.

- **Categorical groups:** Any number of levels in `variable` are
  allowed. Must be a factor variable.

- **Approximate Variance:** Calculations to approximate the variance and
  confidence intervals are available.

- **Allowed Time Values:** Allows both continuous and integer time.

- **Bounded Estimates:** Estimates are not guaranteed to be bounded in
  the 0 to 1 probability range.

- **Monotone Function:** Estimates are not guaranteed to be monotone.

- **Dependencies:** This method relies on the geepack and prodlim
  packages. Additionally requires the eventglm package if
  `censoring_vars` is specified.

Instead of only modeling the outcome mechanism or the treatment
assignment mechanism, both kind of models are required to use this
method. If either of those models are correctly specified, unbiased
estimates will be obtained. In contrast to the
[`"aiptw"`](https://robindenz1.github.io/adjustedCurves/reference/surv_aiptw.md)
method, the "aiptw_pseudo" method uses a generalized estimation equation
(geese) approach to model the outcome mechanism. The model is fit in the
same way as described in the
[`"direct_pseudo"`](https://robindenz1.github.io/adjustedCurves/reference/surv_direct_pseudo.md)
method. Those Direct Standardization based estimates are then
transformed using the previously estimated propensity score. This
results in the doubly-robust property of the method. More information on
this particular method can be found in the original article by Wang
(2018), more information on Pseudo-Values is available in Andersen et
al. (2017) and Andersen and Perme (2010).

When estimating the `geese` model the `ev_time` variable is used as a
factor by default. This results in one coefficient being estimated for
each unique point in time, which can be very slow computationally if
there are a lot of unique points in time and/or the dataset has many
rows. In these cases it is recommended to use `type_time="bs"` or
`type_time="ns"`, which results in the `ev_time` being modeled using
B-Splines or Natural Splines. Simulation studies indicate that there is
little difference in the estimates when an appropriately large number of
`spline_df` is used.

Additionally, covariate-dependent censoring can be accounted for by
using inverse probability of censoring weighted pseudo-values (Binder et
al. 2014) instead of regular pseudo-values (specified using the
`censoring_vars` and `ipcw_method` arguments).

## Value

Adds the following additional objects to the output of the
`adjustedsurv` function:

- `pseudo_values`: The matrix of estimated pseudo-values.

- `geese_model`: The geese model used to make the predictions.

## References

Jixian Wang (2018). "A Simple, Doubly Robust, Efficient Estimator for
Survival Functions Using Pseudo Observations". In: Pharmaceutical
Statistics 17.38-48

James M. Robins and Andrea Rotnitzky (1992). "Recovery of Information
and Adjustment for Dependent Censoring Using Surrogate Markers". In:
AIDS Epidemiology: Methodological Issues. Ed. by Nicholas P. Jewell,
Klaus Dietz, and Vernon T. Farewell. New York: Springer Science +
Business Media, pp. 297-331

Per Kragh Andersen, Elisavet Syriopoulou, and Erik T. Parner (2017).
"Causal Inference in Survival Analysis using Pseudo-Observations". In:
Statistics in Medicine 36, pp. 2669-2681

Per Kragh Andersen and Maja Pohar Perme (2010). "Pseudo-Observations in
Survival Analysis". In: Statistical Methods in Medical Research 19, pp.
71-99

Aris Perperoglou, Willi Sauerbrei, Michal Abrahamowicz, and Matthias
Schmid (2019). "A Review of Spline Function Procedures in R". in: BMC
Medical Research Methodology 19.46, pp. 1-16

Nadine Binder, Thomas A. Gerds, and Per Kragh Andersen (2014). "Pseudo-
Observations for Competing Risks with Covariate Dependent Censoring".
In: Lifetime Data Analysis 20, pp. 303-315

## Author

Jixian Wang supplied the R source code used in the original article,
which was used by Robin Denz to create a generalized version of this
method with additional functionality and improved performance.

## See also

[`geese`](https://rdrr.io/pkg/geepack/man/geese.html),
[`jackknife`](https://rdrr.io/pkg/prodlim/man/jackknife.html),
[`ns`](https://rdrr.io/r/splines/ns.html),
[`bs`](https://rdrr.io/r/splines/bs.html)

## Examples

``` r
library(adjustedCurves)

if (requireNamespace("geepack") & requireNamespace("prodlim")) {

set.seed(42)

# simulate some data as example
sim_dat <- sim_confounded_surv(n=50, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# estimate a treatment assignment model
glm_mod <- glm(group ~ x1 + x3 + x5 + x6, data=sim_dat, family="binomial")

# use it + pseudo values + geese model to calculate adjusted survival curves
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="aiptw_pseudo",
                        outcome_vars=c("x1", "x2", "x3", "x4", "x5", "x6"),
                        treatment_model=glm_mod,
                        conf_int=TRUE,
                        iso_reg=TRUE,
                        force_bounds=TRUE)

# plot the curves
plot(adjsurv, conf_int=TRUE)
}
#> Ignoring unknown labels:
#> • linetype : "Group"
```
