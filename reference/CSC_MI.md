# Cause-Specific Cox Regression with Multiple Imputation

This function can be utilized to perform Cause-Specific Cox Regression
on multiply imputed datasets.

## Usage

``` r
CSC_MI(mids, formula, ...)
```

## Arguments

- mids:

  A `mids` object created using the
  [`mice`](https://amices.org/mice/reference/mice.html) function. This
  replaces the `data` argument in the original function call.

- formula:

  A formula object passed to the
  [`CSC`](https://rdrr.io/pkg/riskRegression/man/CSC.html) function in
  the riskRegression package.

- ...:

  Other arguments which should be passed to the
  [`CSC`](https://rdrr.io/pkg/riskRegression/man/CSC.html) function in
  the riskRegression package.

## Details

A small convenience function to perform CSC regression on multiply
imputed data. It is simply a wrapper around the
[`CSC`](https://rdrr.io/pkg/riskRegression/man/CSC.html) function from
the riskRegression package, because the usual use of `with` is not
supported directly. It returns a `mira` object, which can be passed to
the `outcome_model` argument inside of the
[`adjustedcif`](https://robindenz1.github.io/adjustedCurves/reference/adjustedcif.md)
function when needed. No `pool` method or other functionality is
available.

## Value

A `mira` object containing the CSC regression for every imputed dataset.

## Author

Robin Denz

## See also

[`adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md),
[`CSC`](https://rdrr.io/pkg/riskRegression/man/CSC.html),
[`mice`](https://amices.org/mice/reference/mice.html)

## Examples

``` r
# not run because it would be too slow
# \donttest{
library(adjustedCurves)
library(survival)

if (requireNamespace("riskRegression") & requireNamespace("mice")) {
library(riskRegression)
library(mice)

# simulate some data as example
sim_dat <- sim_confounded_crisk(n=50, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# introduce random missingness in x1 as example
sim_dat$x1 <- ifelse(runif(n=50) < 0.5, sim_dat$x1, NA)

# perform multiple imputation
mids <- mice::mice(data=sim_dat, method="pmm", m=5, printFlag=0)

# use the function
csc_mods <- CSC_MI(mids=mids,
                   formula=Hist(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6 + group
                   )
}
#> Loading required namespace: riskRegression
#> Loading required namespace: mice
#> Registered S3 method overwritten by 'broom':
#>   method        from          
#>   nobs.multinom riskRegression
#> riskRegression version 2025.09.17
#> 
#> Attaching package: ‘mice’
#> The following object is masked from ‘package:stats’:
#> 
#>     filter
#> The following objects are masked from ‘package:base’:
#> 
#>     cbind, rbind
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
# }
```
