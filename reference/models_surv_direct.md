# List of supported models in `surv_direct`

Supported models for the `outcome_model` argument when using
`method="direct"` in the
[`adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md)
function.

## Details

The following models are directly supported in the `outcome_model` in
the
[`surv_direct`](https://robindenz1.github.io/adjustedCurves/reference/surv_direct.md)
function. The first letter in parentheses after the object name is a
group indicator. Below the list there are more information for each
group.

- [`coxph`](https://rdrr.io/pkg/survival/man/coxph.html) \[**A**,
  Required Packages: survival, riskRegression\]

- [`cph`](https://rdrr.io/pkg/rms/man/cph.html) \[**A**, Required
  Packages: rms, survival, riskRegression\]

- [`aalen`](https://rdrr.io/pkg/timereg/man/aalen.html) \[**B**,
  Required Packages: timereg, pec\]

- [`cox.aalen`](https://rdrr.io/pkg/timereg/man/cox.aalen.html) \[**B**,
  Required Packages: timereg, pec\]

- [`selectCox`](https://rdrr.io/pkg/riskRegression/man/selectCox.html)
  \[**B**, Required Packages: riskRegression, pec\]

- [`pecCforest`](https://rdrr.io/pkg/pec/man/pecCforest.html) \[**B**,
  Required Packages: pec\]

- [`pecRpart`](https://rdrr.io/pkg/pec/man/pecRpart.html) \[**B**,
  Required Packages: pec, Bootstrapping not allowed.\]

- [`riskRegression`](https://rdrr.io/pkg/riskRegression/man/riskRegression.html)
  \[**C**, Required Packages: riskRegression\]

- [`prodlim`](https://rdrr.io/pkg/prodlim/man/prodlim.html) \[**C**,
  Required Packages: prodlim, riskRegression\]

- [`psm`](https://rdrr.io/pkg/rms/man/psm.html) \[**C**, Required
  Packages: rms, riskRegression\]

- `flexsurvreg` \[**C**, Required Packages: flexsurv, riskRegression\]

- `flexsurvspline` \[**C**, Required Packages: flexsurv,
  riskRegression\]

- [`ranger`](http://imbs-hl.github.io/ranger/reference/ranger.md)
  \[**C**, Required Packages: ranger, riskRegression\]

- `rfsrc` \[**C**, Required Packages: randomForestSRC, riskRegression\]

- `ARR` \[**C**, Required Packages: riskRegression\]

- [`penalizedS3`](https://rdrr.io/pkg/riskRegression/man/penalizedS3.html)
  \[**C**, Required Packages: penalized, riskRegression\]

- `gbm` \[**C**, Required Packages: gbm, riskRegression\]

- `fit_hal` \[**C**, Required Packages: hal9001, riskRegression\]

- `fitSmoothHazard` \[**C**, Required Packages: casebase,
  riskRegression\]

- [`glm`](https://rdrr.io/r/stats/glm.html) \[**D**, Required Packages:
  stats, pec\]

- [`ols`](https://rdrr.io/pkg/rms/man/ols.html) \[**D**, Required
  Packages: rms, pec\]

- `randomForest` \[**D**, Required Packages: randomForest, pec\]

- `mexhaz` \[**E**, Required Packages: mexhaz\]

- Any model with a fitting S3 prediction method or a valid `predict_fun`
  can be used as well. See below.

**Group A:** The direct adjusted survival probabilities are estimated
directly using the
[`ate`](https://rdrr.io/pkg/riskRegression/man/ate.html) function.
Additional arguments supplied using the `...` syntax are passed to the
[`ate`](https://rdrr.io/pkg/riskRegression/man/ate.html) function. Note
that [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) calls
required to fit the model should be made inside the formula, not saved
elsewhere.  
**Group B:** Predicted survival probabilities are obtained using the
[`predictSurvProb`](https://rdrr.io/pkg/pec/man/predictSurvProb.html)
function. The G-Computation is carried out using those. Additional
arguments supplied using the `...` syntax are passed to the
[`predictSurvProb`](https://rdrr.io/pkg/pec/man/predictSurvProb.html)
function.  
**Group C:** The
[`predictRisk`](https://rdrr.io/pkg/riskRegression/man/predictRisk.html)
function is used to obtain predicted cumulative incidences, which are
then transformed to survival probabilities. Additional arguments
supplied using the `...` syntax are passed to the
[`predictRisk`](https://rdrr.io/pkg/riskRegression/man/predictRisk.html)
function.  
**Group D:** These models are only allowed if there is no censoring.
Predicted survival probabilities are obtained using the `predictProb`
function from the pec package. Additional arguments supplied using the
`...` syntax are passed to the `predictProb` function.  
**Group E:** Custom code is used to obtain predicted survival
probabilities. Additional arguments are not used.

It is sometimes possible to use models even if they are not listed here.
There are two ways to make this work. The first one is to use the models
S3 `predict` method. This works if the `predict` function contains the
arguments `object`, `newdata` and `times` and returns a matrix of
predicted survival probabilities. The matrix should be of size
`nrow(data) * length(times)`, where each row corresponds to a row in the
original dataset and each column to one point in time. The matrix should
contain the survival probabilities predicted by the model given
covariates. If no such `predict` method exists the only option left is
to write your own function which produces the output described above and
supply this function to the `predict_fun` argument.

If you think that some important models are missing from this list,
please file an issue on the official github page with a specific feature
request (URL can be found in the DESCRIPTION file) or contact the
package maintainer directly using the given e-mail address.

## Note

When using outcome models which are not directly supported (either
through the default predict method or a custom `predict_fun`) it might
be necessary to set the `clean_data` argument of the `adjustedsurv`
function to `FALSE`.
