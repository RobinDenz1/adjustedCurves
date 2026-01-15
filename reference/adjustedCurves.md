# Confounder-Adjusted Survival Curves and Cumulative Incidence Functions

***What is this package about?***

This package aims to unite all available adjustments methods for
estimate confounder-adjusted survival curves and cause-specific
confounder-adjusted cumulative incidence functions under one consistent
framework. We try to make the usage of these methods and the calculation
of associated statistics as easy as possible for the user, while still
providing substantial functionality.

***What exactly are adjusted survival curves / adjusted cumulative
incidence functions?***

It is well known that confounding is a serious problem when analyzing
data from non-randomized studies. This is also true when estimating
survival curves or cumulative incidence functions. The aim is to
estimate the population averaged survival probability or cumulative
incidence for some group \\z\\, which would have been observed if every
individual would have been assigned to group \\z\\. For example, the
formal definition for the causal survival curve is:

\$\$S\_{z}(t) = E(I(T_z \> t))\$\$

where \\T_z\\ is the survival time that would have been observed if
treatment \\z\\ was actually administered. See Denz et al. (2023) or Cai
and Van der Laan (2020) for more details

***What features are included in this package?***

This package includes 15 methods to estimate confounder-adjusted
survival curves (single event) and 7 methods to estimate confounder
adjusted cumulative incidence functions (possibly with multiple
competing events). It provides `plot` functions to easily produce highly
customizable and publication-ready graphics. It also allows the user to
easily calculate relevant statistics, such as confidence intervals,
p-values, and adjusted restricted mean survival time estimates. Multiple
Imputation is directly supported.

***What does a typical workflow using this package look like?***

The design of this package is based on the design of the WeightIt
package. It includes two main functions:
[`adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md)
and
[`adjustedcif`](https://robindenz1.github.io/adjustedCurves/reference/adjustedcif.md).
Every implemented adjustment method has their own documentation page
including a small description, code examples, and relevant literature
references. The typical workflow using this package is as follows (1)
estimate confounder-adjusted curves (survival curves or CIFs) using
either `adjustedsurv` or `adjustedcif`, (2) plot those using the S3
`plot` method and sometimes (3) calculate further statistics using
`adjusted_rmst`, `adjusted_rmtl` or `adjusted_curve_test`.

***When should I use `adjustedsurv` and when `adjustedcif`?***

With standard time-to-event data where only one type of event is
possible both the confounder-adjusted survival curves and the
confounder-adjusted cumulative incidence function can be estimated using
the `adjustedsurv` function. While the adjustedsurv function only
estimates the survival, the CIF can simply be calculated by \\1 -
S(t)\\. This transformation from survival curves to CIFs is directly
implemented in the `plot` function (argument `cif`).

When competing risks are present, the cause-specific confounder-adjusted
survival curves can not be estimated in an unbiased way (see for example
Satagopan et al. (2004) for an explanation). The cause-specific
confounder-adjusted cumulative incidence functions however can be
estimated using the `adjustedcif` function.

***What features are missing from this package?***

In a former version, this package included two targeted maximum
likelihood based methods for the estimation of adjusted survival curves
and one for the estimation of adjusted cumulative incidence functions
based on discrete-time data. These methods have been removed because the
survtmle package has been removed from CRAN and there is currently no
other available implementation of these estimators on CRAN. From version
0.10.2 tp 0.11.0 this package then contained an implementation of
targeted maximum likelihood estimation for continuous time-to-event data
(a wrapper function for the concrete package). This is currently
unavailable due to the package being removed from CRAN, but will
probably be supported again soon.

This package also currently does not support time-varying treatments or
covariates. It also does not support left-censoring, interval-censoring
or left-truncation. These features may be added in future releases.

***What if the variable of interest is continuous?***

All methods in this package are designed strictly for categorical
variables of interest. If the variable of interest is continuous the
user has to manually categorize the variable and save it as a factor
first. This is, however, generally discouraged because artificial
categorization may lead to bias or misleading results. To face this
issue we have developed the `contsurvplot` R package which implements
multiple plotting methods to visualize the (causal) effect of a
continuous variable when analyzing survival data (Denz & Timmesfeld
2023).

***Where can I get more information?***

The documentation pages contain a lot of information, relevant examples
and literature references. Additional examples can be found in the
vignette of this package, which can be accessed using
[`vignette(topic="adjustedCurves", package="adjustedCurves")`](https://robindenz1.github.io/adjustedCurves/articles/adjustedCurves.md)
or in the main paper associated with this package (Denz et al. 2023). We
are also working on a separate article on this package that is going to
be published in a peer-reviewed journal.

***I want to suggest a new feature / I want to report a bug. Where can I
do this?***

Bug reports, suggestions and feature requests are highly welcome. Please
file an issue on the official github page or contact the author directly
using the supplied e-mail address.

## References

Robin Denz, Renate Klaaßen-Mielke, and Nina Timmesfeld (2023). "A
Comparison of Different Methods to Adjust Survival Curves for
Confounders". In: Statistics in Medicine 42.10, pp. 1461-1479.

Robin Denz, Nina Timmesfeld (2023). "Visualizing the (Causal) Effect of
a Continuous Variable on a Time-To-Event Outcome". In: Epidemiology
34.5, pp. 652-660.

Weixin Cai and Mark J. van der Laan (2020). "One-Step Targeted Maximum
Likelihood Estimation for Time-To-Event Outcomes". In: Biometrics 76,
pp. 722-733

J. M. Satagopan, L. Ben-Porat, M. Berwick, M. Robson, D. Kutler, and A.
D. Auerbach (2004). "A Note on Competing Risks in Survival Data
Analysis". In: British Journal of Cancer 91, pp. 1229-1235.

## Author

Robin Denz, \<robin.denz@rub.de\>
