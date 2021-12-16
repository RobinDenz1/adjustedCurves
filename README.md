<!-- badges: start -->
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![](https://www.r-pkg.org/badges/version/badger?color=green)](https://cran.r-project.org/package=badger)
[![](http://cranlogs.r-pkg.org/badges/grand-total/adjustedCurves?color=green)](https://cran.r-project.org/package=adjustedCurves)
[![R-CMD-check](https://github.com/RobinDenz1/adjustedCurves/workflows/R-CMD-check/badge.svg)](https://github.com/RobinDenz1/adjustedCurves/actions)
<!-- badges: end -->

# adjustedCurves

Author: Robin Denz

## Description

`adjustedCurves` is an R-Package which can be used to calculate and plot confounder-adjusted survival curves + confidence intervals as well as cause-specific confounder-adjusted cumulative incidence functions + confidence intervals using a variety of methods.
It provides an convenient wrapper around existing R-Packages on the topic and adds additional methods and functionality on top of it.
Those additional features include the calculation of adjusted restricted mean survival time and testing the equality of two confounder-adjusted survival curves.

Detailed descriptions of each method can be found in the literature cited in the documentation. 

## Installation

Currently this package is not available on CRAN, but can be installed easily using the `devtools` R-Package:

```
library(devtools)

devtools::install_github("https://github.com/RobinDenz1/adjustedCurves")
```

or the `remotes` R-Package:

```
library(remotes)

remotes::install_github("https://github.com/RobinDenz1/adjustedCurves")
```

## Issues

If you encounter any bugs or have any specific feature requests, please file an issue.

## Examples

This minimal example shows how to calculate adjusted survival curves using *Direct Adjustment* with this package:

```
library(adjustedCurves)

# simulate some data as example
sim_dat <- sim_confounded_surv(n=500)
sim_dat$group <- as.factor(sim_dat$group)

# take a look at the data
head(sim_dat)

# estimate a cox-regression for the outcome
cox_mod <- coxph(Surv(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6 + group,
                 data=sim_dat, x=T)


# use it to calculate adjusted survival curves
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="direct",
                        outcome_model=cox_mod,
                        conf_int=T)

# plot the curves
plot(adjsurv, draw_ci=F)

# also plot the confidence intervals
plot(adjsurv)
```
Here is an example of how to calculate adjusted survival curves using *Inverse Probability of Treatment Weighting*:
```
# estimate a treatment assignment model
glm_mod <- glm(group ~ x1 + x2 + x3 + x4 + x5 + x6, data=sim_dat,
               family="binomial"(link="logit"))


# use it to calculate adjusted survival curves
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="iptw_km",
                        treatment_model=glm_mod,
                        conf_int=T)

# plot the curves
plot(adjsurv, draw_ci=F)
```
To test whether the two adjusted survival curves are equal, the `adjustedsurv` call has to be made with `bootstrap=TRUE`:
```
# use it to calculate adjusted survival curves
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="direct",
                        outcome_model=cox_mod,
                        conf_int=T,
                        bootstrap=T,
                        n_boot=1000)
                        
 adj_test <- test_curve_equality(adjsurv, from=0, to=1.2)
 print(adj_test)
```

More examples can be found in the documentation and the vignettes.

## Citation
Please cite this R-Package using:
MY PAPER

You should also cite the paper describing the method you used. The respective literature can be found in the documentation.

## License

© 2021-2021 Robin Denz

The contents of this repository are distributed under the GNU General Public License. You can find the full text of this License in this github repository. Alternatively, see <http://www.gnu.org/licenses/>.


