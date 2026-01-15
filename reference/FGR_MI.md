# Fine & Gray Model with Multiple Imputation

This function can be utilized to calculate Fine & Gray models for
multiply imputed datasets.

## Usage

``` r
FGR_MI(mids, formula, cause=1, ...)
```

## Arguments

- mids:

  A `mids` object created using the
  [`mice`](https://amices.org/mice/reference/mice.html) function. This
  replaces the `data` argument in the original function call.

- formula:

  A formula object passed to the
  [`FGR`](https://rdrr.io/pkg/riskRegression/man/FGR.html) function in
  the riskRegression package.

- cause:

  The failure type of interest. Defaults to 1.

- ...:

  Other arguments which should be passed to the
  [`FGR`](https://rdrr.io/pkg/riskRegression/man/FGR.html) function in
  the riskRegression package.

## Details

A small convenience function to calculate Fine & Gray models for
multiply imputed data. It is simply a wrapper around the
[`FGR`](https://rdrr.io/pkg/riskRegression/man/FGR.html) function from
the riskRegression package, because the usual use of `with` is not
supported directly. It returns a `mira` object, which can be passed to
the `outcome_model` argument inside of the
[`adjustedcif`](https://robindenz1.github.io/adjustedCurves/reference/adjustedcif.md)
function when needed. No `pool` method or other functionality is
available.

## Value

A `mira` object containing the FGR regression for every imputed dataset.

## Author

Robin Denz

## See also

[`adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md)

## Examples

``` r
# not run because it would be too slow
# \donttest{
library(adjustedCurves)
library(survival)

if (requireNamespace("riskRegression") & requireNamespace("prodlim") &
    requireNamespace("mice")) {

library(riskRegression)
library(mice)
library(prodlim)

# simulate some data as example
sim_dat <- sim_confounded_crisk(n=50, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# introduce random missingness in x1 as example
sim_dat$x1 <- ifelse(runif(n=50) < 0.5, sim_dat$x1, NA)

# perform multiple imputation
mids <- mice::mice(data=sim_dat, method="pmm", m=5, printFlag=FALSE)

# use the function
fgr_mods <- FGR_MI(mids=mids,
                   formula=Hist(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6 + group,
                   cause=1)
}
# }
```
