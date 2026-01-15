# Direct Adjusted Survival Curves

This page explains the details of estimating confounder-adjusted
survival curves using a previously fit Cox-Regression model for single
event time-to-event data using Direct Standardization (`method="direct"`
in the
[`adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md)
function). All regular arguments of the `adjustedsurv` function can be
used. Additionally, the `outcome_model` argument has to be specified in
the `adjustedsurv` call. Further arguments specific to this method are
listed below.

## Arguments

- outcome_model:

  \[**required**\] Must be a previously fit model object including
  `variable` as independent variable. Apart from the classic `coxph`
  model this function also supports a variety of other models. See
  [`models_surv_direct`](https://robindenz1.github.io/adjustedCurves/reference/models_surv_direct.md)
  for a list of supported model objects and some more details.

- verbose:

  Whether to print estimation information of the `ate` function in the
  riskRegression package. Ignored if `outcome_model` is not a `coxph`
  object. Defaults to `FALSE`.

- predict_fun:

  A function which should be used to calculate the predicted survival
  probabilities given covariates and some points in time. This argument
  only needs to be specified if the kind of model supplied in the
  `outcome_model` is not directly supported. See
  [`models_surv_direct`](https://robindenz1.github.io/adjustedCurves/reference/models_surv_direct.md)
  for more information. Defaults to `NULL`.

- ...:

  Further arguments passed to `ate` if `outcome_model` is a `coxph`
  object. Otherwise the additional arguments are passed to the
  respective `predict` method. See
  [`models_surv_direct`](https://robindenz1.github.io/adjustedCurves/reference/models_surv_direct.md)
  for more information.

## Details

- **Type of Adjustment:** Requires a model describing the outcome
  mechanism. See
  [`models_surv_direct`](https://robindenz1.github.io/adjustedCurves/reference/models_surv_direct.md)
  for a list of supported model objects and some more details.

- **Doubly-Robust:** Estimates are not Doubly-Robust.

- **Categorical groups:** Any number of levels in `variable` are
  allowed. Must be a factor variable.

- **Approximate Variance:** Calculations to approximate the variance and
  confidence intervals are available only if `outcome_model` is a
  `coxph` object. The
  [`ate`](https://rdrr.io/pkg/riskRegression/man/ate.html) function is
  used for the calculation in that case. Bootstrap confidence intervals
  can however be calculated with all supported models. See
  [`?adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md)
  for more information on bootstrapping.

- **Allowed Time Values:** Allows both continuous and integer time.

- **Bounded Estimates:** Estimates are guaranteed to be bounded in the 0
  to 1 probability range.

- **Monotone Function:** Estimates are guaranteed to be monotone.

- **Dependencies:** This method relies on the riskRegression package.
  Depending on `outcome_model` other packages might be needed. See
  [`models_surv_direct`](https://robindenz1.github.io/adjustedCurves/reference/models_surv_direct.md)
  for more details.

This method works by executing the following steps: (1) First a model is
fitted which describes the outcome mechanism (time-to-event). Next (2)
multiple copies of the original dataset are created, one for each
possible level of the `variable` of interest. (3) The `variable` is then
set to one level for all observations in each dataset. (4) The model is
used to predict the survival probabilities at some points in time T for
each observation in all dataset copies. (5) Those estimated
probabilities are averaged for each dataset at each point in time,
resulting in adjusted survival probabilities for all levels of the group
variable at the specified points in time.

In the literature this method is sometimes called "Direct
Standardization", "Corrected Group-Prognosis", "G-Computation" or
"G-Formula". If the model in step (1) is "correct"" this method will
produce unbiased estimates of the counterfactual survival curves. A
model can be called a "correct" model in this context if it can be used
to produce unbiased estimates of the true (but unknown) individual
survival probabilities given covariates. When used properly this is one
of the most efficient methods. More information can be found in the
literature listed in the references. The most popular model for
describing the outcome mechanism in a time-to-event context is the
Cox-regression model
([`coxph`](https://rdrr.io/pkg/survival/man/coxph.html)). This function
however also supports a variety of other models.

## Value

Adds the following additional objects to the output of the
`adjustedsurv` function:

- `ate_object`: The object returned by the `ate` function.

## References

I-Ming Chang, Rebecca Gelman, and Marcello Pagano (1982). "Corrected
Group Prognostic Curves and Summary Statistics". In: Journal of Chronic
Diseases 35, pp. 669-674

Robert W. Makuch (1982). "Adjusted Survival Curve Estimation Using
Covariates". In: Journal of Chronic Diseases 35.6, pp. 437-443

Xu Zhang, Fausto R. Loberiza, John P. Klein, and Mei-Jie Zhang (2007).
"A SAS Macro for Estimation of Direct Adjusted Survival Curves Based on
a Stratified Cox Regression Model". In: Computer Methods and Programs in
Biomedicine 88, pp. 95-101

## Author

The function itself was written by Robin Denz. When using `coxph` models
however, this function is just a wrapper around the
[`ate`](https://rdrr.io/pkg/riskRegression/man/ate.html) function, which
was written by other people. See
[`?ate`](https://rdrr.io/pkg/riskRegression/man/ate.html) for more
information.

## See also

[`models_surv_direct`](https://robindenz1.github.io/adjustedCurves/reference/models_surv_direct.md),
[`ate`](https://rdrr.io/pkg/riskRegression/man/ate.html),
[`coxph`](https://rdrr.io/pkg/survival/man/coxph.html)

## Examples

``` r
library(adjustedCurves)
library(survival)

if (requireNamespace("riskRegression")) {

library(riskRegression)

set.seed(42)

# simulate some data as example
sim_dat <- sim_confounded_surv(n=50, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# estimate a cox-regression for the outcome
# NOTE: some authors also recommend using strata(group), see FAQ vignette
cox_mod <- coxph(Surv(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6 + group,
                 data=sim_dat, x=TRUE)

# use it to calculate adjusted survival curves
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="direct",
                        outcome_model=cox_mod,
                        conf_int=FALSE)

# plot the curves
plot(adjsurv)

# not run to avoid dependency on flexsurv and mice too slow
if (interactive()) {
## using a flexsurv() model, this requires the 'fleysurv' package
mod_flexsurvreg <- flexsurvreg(Surv(time, event) ~ group + x1 + x2 + x5 + x6,
                               data=sim_dat, dist="gengamma")

# using it to calculate the adjusted survival curves
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="direct",
                        outcome_model=mod_flexsurvreg,
                        conf_int=FALSE)

# plot using steps=FALSE to draw them as smooth functions, since
# they were estimated using a parametric model
plot(adjsurv, steps=FALSE)
}

# \donttest{
## using multiple imputation
if (requireNamespace("mice")) {
library(mice)

# introduce random missingness in x1 as example
# NOTE: This is only done as an example, in reality you would
#       already have missing data, not introduce it yourself.
sim_dat$x1 <- ifelse(runif(n=50) < 0.5, sim_dat$x1, NA)

# perform multiple imputation
mids <- mice::mice(data=sim_dat, method="pmm", m=5, printFlag=FALSE)

# fit model for each imputed dataset
mira <- with(mids, coxph(Surv(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6 + group,
                         x=TRUE))

# calculate adjusted survival curves on imputed data
adj <- adjustedsurv(data=mids,
                    variable="group",
                    ev_time="time",
                    event="event",
                    method="direct",
                    outcome_model=mira)
plot(adj)
}
# }
}
#> Ignoring unknown labels:
#> • linetype : "Group"
#> • fill : "Group"
```
