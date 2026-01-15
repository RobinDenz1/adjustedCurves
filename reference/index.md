# Package index

## Adjusted survival curves

Estimate and plot confounder adjusted survival curves using various
methods

- [`adjustedsurv()`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md)
  : Estimate Confounder-Adjusted Survival Curves

- [`plot(`*`<adjustedsurv>`*`)`](https://robindenz1.github.io/adjustedCurves/reference/plot.adjustedsurv.md)
  : Plot Confounder-Adjusted Survival Curves

- [`surv_direct`](https://robindenz1.github.io/adjustedCurves/reference/surv_direct.md)
  : Direct Adjusted Survival Curves

- [`models_surv_direct`](https://robindenz1.github.io/adjustedCurves/reference/models_surv_direct.md)
  :

  List of supported models in `surv_direct`

- [`surv_direct_pseudo`](https://robindenz1.github.io/adjustedCurves/reference/surv_direct_pseudo.md)
  : Direct Adjusted Survival Curves using Pseudo-Values

- [`surv_iptw_km`](https://robindenz1.github.io/adjustedCurves/reference/surv_iptw_km.md)
  : Inverse Probability of Treatment Weighted Kaplan-Meier estimates

- [`surv_iptw_cox`](https://robindenz1.github.io/adjustedCurves/reference/surv_iptw_cox.md)
  : Inverse Probability of Treatment Weighted Survival using
  Cox-Regression

- [`surv_iptw_pseudo`](https://robindenz1.github.io/adjustedCurves/reference/surv_iptw_pseudo.md)
  : Inverse Probability of Treatment Weighted Survival Estimates using
  Pseudo-Values

- [`surv_matching`](https://robindenz1.github.io/adjustedCurves/reference/surv_matching.md)
  : Using Propensity-Score Matching to Calculate Adjusted Survival
  Curves

- [`surv_emp_lik`](https://robindenz1.github.io/adjustedCurves/reference/surv_emp_lik.md)
  : Empirical Likelihood Estimation Survival Curves

- [`surv_aiptw`](https://robindenz1.github.io/adjustedCurves/reference/surv_aiptw.md)
  : Augmented Inverse Probability of Treatment Weighted Survival Curves

- [`surv_aiptw_pseudo`](https://robindenz1.github.io/adjustedCurves/reference/surv_aiptw_pseudo.md)
  : Augmented Inverse Probability of Treatment Weighted Survival Curves
  using Pseudo-Values

- [`surv_strat_amato`](https://robindenz1.github.io/adjustedCurves/reference/surv_strat_amato.md)
  : Adjusted Survival Curves for Categorical Confounders using the
  Method by Amato (1988)

- [`surv_strat_nieto`](https://robindenz1.github.io/adjustedCurves/reference/surv_strat_nieto.md)
  : Adjusted Survival Curves for Categorical Confounders using the
  Method by Gregory (1988) and Nieto & Coresh (1996)

- [`surv_strat_cupples`](https://robindenz1.github.io/adjustedCurves/reference/surv_strat_cupples.md)
  : Adjusted Survival Curves for Categorical Confounders using the
  Method by Cupples et al. (1995)

- [`surv_iv_2SRIF`](https://robindenz1.github.io/adjustedCurves/reference/surv_iv_2SRIF.md)
  : Instrumental Variable based Survival Curve Estimation using the Two
  Stage Residual Inclusion method with a Frailty Term (2SRI-F)

- [`surv_prox_iptw`](https://robindenz1.github.io/adjustedCurves/reference/surv_prox_iptw.md)
  : Proximal Inverse Probability of Treatment Weighted Survival Curve
  Estimates

- [`surv_prox_aiptw`](https://robindenz1.github.io/adjustedCurves/reference/surv_prox_aiptw.md)
  : Proximal Augmented Inverse Probability of Treatment Weighted
  Survival Curve Estimates

- [`surv_km`](https://robindenz1.github.io/adjustedCurves/reference/surv_km.md)
  : Group-Specific Kaplan-Meier Survival Curves

## Adjusted cumulative incidence functions

Estimate and plot confounder adjusted cumulative incidence functions
from competing events data using various methods

- [`adjustedcif()`](https://robindenz1.github.io/adjustedCurves/reference/adjustedcif.md)
  : Estimate Cause-Specific Confounder-Adjusted Cumulative Incidence
  Functions

- [`plot(`*`<adjustedcif>`*`)`](https://robindenz1.github.io/adjustedCurves/reference/plot.adjustedcif.md)
  : Plot Confounder-Adjusted Cumulative Incidence Functions

- [`cif_direct`](https://robindenz1.github.io/adjustedCurves/reference/cif_direct.md)
  : Direct Adjusted Cumulative Incidence Functions

- [`models_cif_direct`](https://robindenz1.github.io/adjustedCurves/reference/models_cif_direct.md)
  :

  List of supported models in `cif_direct`

- [`cif_direct_pseudo`](https://robindenz1.github.io/adjustedCurves/reference/cif_direct_pseudo.md)
  : Direct Adjusted CIFs using Pseudo-Values

- [`cif_iptw`](https://robindenz1.github.io/adjustedCurves/reference/cif_iptw.md)
  : Inverse Probability of Treatment Weighted CIFs

- [`cif_iptw_pseudo`](https://robindenz1.github.io/adjustedCurves/reference/cif_iptw_pseudo.md)
  : Inverse Probability of Treatment Weighted CIFs using Pseudo-Values

- [`cif_matching`](https://robindenz1.github.io/adjustedCurves/reference/cif_matching.md)
  : Using Propensity-Score Matching to Calculate Adjusted CIFs

- [`cif_aiptw`](https://robindenz1.github.io/adjustedCurves/reference/cif_aiptw.md)
  : Augmented Inverse Probability of Treatment Weighted CIFs

- [`cif_aiptw_pseudo`](https://robindenz1.github.io/adjustedCurves/reference/cif_aiptw_pseudo.md)
  : Augmented Inverse Probability of Treatment Weighted CIFs using
  Pseudo-Values

- [`cif_aalen_johansen`](https://robindenz1.github.io/adjustedCurves/reference/cif_aalen_johansen.md)
  : Group-Specific Aalen-Johansen CIFs

## Statistics and Group Comparisons

Various Statistics and ways to compare groups based on adjusted survival
curves or CIFs

- [`adjusted_curve_diff()`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_curve_diff.md)
  [`adjusted_curve_ratio()`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_curve_diff.md)
  : Estimate the difference between or the ratio of two
  Confounder-Adjusted Survival Curves or CIFs

- [`adjusted_rmst()`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_rmst.md)
  : Estimate Confounder-Adjusted Restricted Mean Survival Times

- [`adjusted_rmtl()`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_rmtl.md)
  : Estimate Confounder-Adjusted Restricted Mean Time Lost

- [`adjusted_surv_quantile()`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_surv_quantile.md)
  : Estimate Confounder-Adjusted Survival Time Quantiles

- [`adjusted_curve_test()`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_curve_test.md)
  : Test if there is a difference between two Confounder-Adjusted
  Survival Curves or CIFs

- [`print(`*`<curve_test>`*`)`](https://robindenz1.github.io/adjustedCurves/reference/print.curve_test.md)
  :

  Print Method for `curve_test` Objects

- [`plot(`*`<curve_test>`*`)`](https://robindenz1.github.io/adjustedCurves/reference/plot.curve_test.md)
  :

  Plot Method for `curve_test` Objects

- [`plot_curve_diff()`](https://robindenz1.github.io/adjustedCurves/reference/plot_curve_diff.md)
  [`plot_curve_ratio()`](https://robindenz1.github.io/adjustedCurves/reference/plot_curve_diff.md)
  : Plot the Difference Between or the Ratio of Two Adjusted Survival
  Curves or CIFs

- [`plot_rmst_curve()`](https://robindenz1.github.io/adjustedCurves/reference/plot_rmst_curve.md)
  : Plot Adjusted Restricted Mean Survival Time Curves

- [`plot_rmtl_curve()`](https://robindenz1.github.io/adjustedCurves/reference/plot_rmtl_curve.md)
  : Plot Adjusted Restricted Mean Time Lost Curves

## Misc

- [`adjustedCurves-package`](https://robindenz1.github.io/adjustedCurves/reference/adjustedCurves.md)
  : Confounder-Adjusted Survival Curves and Cumulative Incidence
  Functions

- [`as_ggsurvplot_df()`](https://robindenz1.github.io/adjustedCurves/reference/as_ggsurvplot_df.md)
  :

  Extract a `data.frame` containing the estimated survival curves from a
  `adjustedsurv` object

- [`sim_confounded_surv()`](https://robindenz1.github.io/adjustedCurves/reference/sim_confounded_surv.md)
  : Simulate Survival Data with Confounders

- [`sim_confounded_crisk()`](https://robindenz1.github.io/adjustedCurves/reference/sim_confounded_crisk.md)
  : Simulate Competing Risks Data with Confounders

- [`CSC_MI()`](https://robindenz1.github.io/adjustedCurves/reference/CSC_MI.md)
  : Cause-Specific Cox Regression with Multiple Imputation

- [`FGR_MI()`](https://robindenz1.github.io/adjustedCurves/reference/FGR_MI.md)
  : Fine & Gray Model with Multiple Imputation
