# adjustedCurves

Author: Robin Denz

## Description

`adjustedCurves` is an R-Package which can be used to calculate and plot confounder-adjusted survival curves as well as their confidence intervals using a variety of methods.
It provides an convenient wrapper around existing R-Packages on the topic and adds additional methods and functionality on top of it.
Those additional features include the calculation of adjusted restricted mean survival time and testing the equality of two confounder-adjusted survival curves.

## Installation

Currently this package is not available on Cran, but can be installed easily using the following code:

```
library(devtool)

install_github("https://github.com/RobinDenz1/adjustedCurves")
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
                        sd=T)

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
                        sd=T)

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
                        sd=T,
                        bootstrap=T,
                        n_boot=1000)
                        
 adj_test <- adjustedsurv_test(adjsurv, from=0, to=1.2)
 print(adj_test)
```



