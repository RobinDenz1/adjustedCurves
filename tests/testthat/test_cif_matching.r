library(Matching)

set.seed(42)

sim_dat <- adjustedCurves::sim_confounded_surv(n=50)
sim_dat$event[sim_dat$event==1] <- sample(c(1, 2), size=sum(sim_dat$event),
                                          replace=TRUE)

# estimate propensity score
mod <- glm(group ~ x1 + x2 + x3 + x4 + x5 + x6, data=sim_dat,
           family="binomial")
ps_score <- mod$fitted.values

## Just check if function throws any errors
test_that("2 treatments, no conf_int, no boot", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="matching",
                                           conf_int=FALSE,
                                           treatment_model=mod,
                                           cause=1), NA)
})

test_that("2 treatments, no conf_int, no boot, with ps_score", {
  expect_error(adjustedCurves::adjustedcif(data=sim_dat,
                                           variable="group",
                                           ev_time="time",
                                           event="event",
                                           method="matching",
                                           conf_int=FALSE,
                                           bootstrap=FALSE,
                                           n_boot=2,
                                           treatment_model=ps_score,
                                           cause=1), NA)
})
