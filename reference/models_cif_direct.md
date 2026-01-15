# List of supported models in `cif_direct`

Supported models for the `outcome_model` argument when using
`method="direct"` in the
[`adjustedcif`](https://robindenz1.github.io/adjustedCurves/reference/adjustedcif.md)
function.

## Details

The following models are directly supported in the `outcome_model` in
the
[`cif_direct`](https://robindenz1.github.io/adjustedCurves/reference/cif_direct.md)
function. The first letter in parentheses after the object name is a
group indicator. Below the list there are more information for each
group.

- [`CSC`](https://rdrr.io/pkg/riskRegression/man/CSC.html) \[**A**,
  Required Packages: riskRegression\]

- [`FGR`](https://rdrr.io/pkg/riskRegression/man/FGR.html) \[**B**,
  Required Packages: riskRegression\]

- [`riskRegression`](https://rdrr.io/pkg/riskRegression/man/riskRegression.html)
  \[**B**, Required Packages: riskRegression\]

- [`prodlim`](https://rdrr.io/pkg/prodlim/man/prodlim.html) \[**B**,
  Required Packages: prodlim, riskRegression\]

- `rfsrc` \[**B**, Required Packages: randomForestSRC, riskRegression\]

- `ARR` \[**B**, Required Packages: riskRegression\]

- `fit_hal` \[**B**, Required Packages: hal9001, riskRegression\]

- `fastCrr` \[**C**, Required Packages: fastcmprsk\]

- [`comp.risk`](https://rdrr.io/pkg/timereg/man/comp.risk.html) \[**C**,
  Required Packages: timereg\]

- Any model with a fitting S3 prediction method or a valid `predict_fun`
  can be used as well. See below.

**Group A:** The direct adjusted cumulative incidences are estimated
directly using the
[`ate`](https://rdrr.io/pkg/riskRegression/man/ate.html) function.
Additional arguments supplied using the `...` syntax are passed to the
[`ate`](https://rdrr.io/pkg/riskRegression/man/ate.html) function.  
**Group B:** The
[`predictRisk`](https://rdrr.io/pkg/riskRegression/man/predictRisk.html)
function is used to obtain predicted cumulative incidences, which are
then used in the G-Computation step. Additional arguments supplied using
the `...` syntax are passed to the
[`predictRisk`](https://rdrr.io/pkg/riskRegression/man/predictRisk.html)
function.  
**Group C:** Custom code is used to do the estimation. Additional
arguments supplied using the `...` syntax are currently not supported.  

It is sometimes possible to use models even if they are not listed here.
There are two ways to make this work. The first one is to use the models
S3 `predict` method. This works if the `predict` function contains the
arguments `object`, `newdata`, `times` and `cause` and returns a matrix
of predicted cause-specific cumulative incidences. The matrix should be
of size `nrow(data) * length(times)`, where each row corresponds to a
row in the original dataset and each column to one point in time. The
matrix should contain the cause-specific cumulative incidences predicted
by the model given covariates. If no such `predict` method exists the
only option left is to write your own function which produces the output
described above and supply this function to the `predict_fun` argument.

If you think that some important models are missing from this list,
please file an issue on the official github page with a specific feature
request (URL can be found in the DESCRIPTION file) or contact the
package maintainer directly using the given e-mail address.

## Note

When using outcome models which are not directly supported (either
through the default predict method or a custom `predict_fun`) it might
be necessary to set the `clean_data` argument of the `adjustedcif`
function to `FALSE`.
