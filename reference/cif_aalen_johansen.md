# Group-Specific Aalen-Johansen CIFs

This page explains the details of estimating standard Aalen-Johansen
cumulative incidence functions, stratified by the group variable
(`method="aalen_johansen"` in the
[`adjustedcif`](https://robindenz1.github.io/adjustedCurves/reference/adjustedcif.md)
function). All regular arguments of the `adjustedcif` function can be
used. Further arguments specific to this method are listed below.

NO adjustment for any confounders are made. This function is included
only for reference and should not be used when confounder adjusted CIFs
are desired.

## Arguments

- ...:

  Further arguments passed to
  [`cuminc`](https://rdrr.io/pkg/cmprsk/man/cuminc.html).

## Details

- **Type of Adjustment:** NO adjustments are made. This is just a
  stratified Aalen-Johansen estimator.

- **Doubly-Robust:** Estimates are not Doubly-Robust.

- **Categorical groups:** Any number of levels in `variable` are
  allowed. Must be a factor variable.

- **Approximate Variance:** Calculations to approximate the variance and
  confidence intervals are available.

- **Allowed Time Values:** Allows both continuous and integer time.

- **Bounded Estimates:** Estimates are guaranteed to be bounded in the 0
  to 1 probability range.

- **Monotone Function:** Estimates are guaranteed to be monotone.

- **Dependencies:** This method relies on the the cmprsk package.

This function is just a convenient wrapper around the
[`cuminc`](https://rdrr.io/pkg/cmprsk/man/cuminc.html) function. See
[`?cuminc`](https://rdrr.io/pkg/cmprsk/man/cuminc.html) or the cited
literature for more details.

## Value

Adds the following additional objects to the output of the
`adjustedsurv` function:

- `cuminc_object`: The object returned by the `cuminc` function.

## References

Odd O. Aalen and Søren Johansen (1978). "An Empirical Transition Matrix
for Non-Homogeneous Markov Chains Based on Censored Observations". In:
Scandinavian Journal of Statistics 5.3, pp. 141-150

## Author

The wrapper function was written by Robin Denz, the `cuminc` function
(which this wrapper is build around) was written by other people. See
[`?cuminc`](https://rdrr.io/pkg/cmprsk/man/cuminc.html) for more
details.

## See also

[`cuminc`](https://rdrr.io/pkg/cmprsk/man/cuminc.html)

## Examples

``` r
library(adjustedCurves)

if (requireNamespace("cmprsk")) {

library(cmprsk)

set.seed(42)

# simulate some data as example
sim_dat <- sim_confounded_crisk(n=50, max_t=5)
sim_dat$group <- as.factor(sim_dat$group)

# calculate un-adjusted aalen-johansen estimates
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      cause=1,
                      method="aalen_johansen")

# plot the curves
plot(adjcif)
}
#> Ignoring unknown labels:
#> • linetype : "Group"
#> • fill : "Group"
```
