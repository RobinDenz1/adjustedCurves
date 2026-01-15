# Extract a `data.frame` containing the estimated survival curves from a `adjustedsurv` object

A small convenience function to extract the most important quantities
from an `adjustedsurv` object. The resulting `data.frame` is structured
according to the format required by the `ggsurvplot_df` function of the
survminer package, making it easy to use the `ggsurvplot_df` function.

## Usage

``` r
as_ggsurvplot_df(adjsurv)
```

## Arguments

- adjsurv:

  An object of class `adjustedsurv` created by the
  [`adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md)
  function.

## Value

Returns a `data.frame` containing the required information, extracted
from the `adjustedsurv` object.

## Author

Robin Denz

## See also

[`adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md),
[`plot.adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/plot.adjustedsurv.md)

## Examples

``` r
library(adjustedCurves)
library(survival)

set.seed(42)

# simulate some example data
sim_dat <- sim_confounded_surv(n=50, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# treatment assignment model
glm_mod <- glm(group ~ x2 + x3 + x5 + x6, data=sim_dat, family="binomial")

# estimate some adjusted survival curves
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="iptw_km",
                        treatment_model=glm_mod,
                        conf_int=TRUE,
                        bootstrap=FALSE)

# extract info
df <- as_ggsurvplot_df(adjsurv)

# not run here to avoid dependency on survminer
if (interactive()) {
# plot using survminer, requires the 'survminer' package
ggsurvplot_df(df)
}
```
