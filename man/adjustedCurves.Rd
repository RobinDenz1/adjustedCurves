\name{adjustedCurves-package}
\alias{adjustedCurves-package}

\title{
Confounder-Adjusted Survival Curves and Cumulative Incidence Functions
}
\description{

\strong{\emph{What is this package about?}}

This package aims to unite all available adjustments methods for calculating confounder-adjusted survival curves and confounder-adjusted cumulative incidence functions under one consistent framework. We try to make the usage of these methods and the calculation of associated statistics as easy as possible for the user, while still providing substantial functionality.

\strong{\emph{What exactly are adjusted survival curves / adjusted cumulative incidence functions?}}

It is well known that confounding is a serious problem when analyzing data from non-randomized studies. This is also true when calculating survival curves or cumulative incidence functions. The aim is to estimate the population averaged survival probability or cumulative incidence for some group \eqn{z}, which would have been observed if every individual would have been assigned to group \eqn{z}. A formal definition of the estimand can be found in Cai & van der Laan (2020) on page 723.

\strong{\emph{What features are included in this package?}}

This package includes 14 methods to calculate confounder-adjusted survival curves (single event) and eight methods to calculate confounder adjusted cumulative incidence functions (possibly with multiple competing events). It provides \code{plot} functions to easily produce highly customizable and publication-ready graphics. It also allows the user to easily calculate relevant statistics, such as confidence intervals, p-values, and adjusted restricted mean survival time estimates. Multiple Imputation is directly supported.

\strong{\emph{What does a typical workflow using this package look like?}}

The design of this package is based on the design of the \pkg{WeightIt} package. It includes two main functions: \code{\link[adjustedCurves]{adjustedsurv}} and \code{\link[adjustedCurves]{adjustedcif}}. Every implemented adjustment method has their own documentation page including a small description, code examples, and relevant literature references. The typical workflow using this package is as follows (1) calculate confounder-adjusted curves (survival curves or CIFs) using either \code{adjustedsurv} or \code{adjustedcif}, (2) plot those using the S3 \code{plot} method and sometimes (3) calculate further statistics using \code{adjusted_rmst} or \code{test_curve_equality}.

\strong{\emph{When should I use \code{adjustedsurv} and when \code{adjustedcif}?}}

With standard time-to-event data where only one type of event is possible both the confounder-adjusted survival curves and the confounder-adjusted cumulative incidence function can be estimated using the \code{adjustedsurv} function. While the adjustedsurv function only estimates the survival, the CIF can simply be calculated by \eqn{1 - S(t)}. This transformation from survival curves to CIFs is directly implemented in the \code{plot} function (argument \code{cif}).

When competing risks are present, the cause-specific confounder-adjusted survival curves can not be estimated in an unbiased way (see for example Satagopan et al. (2004) for an explanation). The cause-specific confounder-adjusted cumulative incidence functions however can be estimated using the \code{adjustedcif} function.

\strong{\emph{I want to suggest a new feature / I want to report a bug. Where can I do this?}}

Bug reports, suggestions and feature requests are highly welcome. Please file an issue on the official github page or contact the author directly using the supplied e-mail address.
}
\references{
Weixin Cai & Mark J. van der Laan. "One-Step Targeted Maximum Likelihood Estimation for Time-To-Event Outcomes". In: Biometrics 76 (2020), pp. 722-733.

J. M. Satagopan; L. Ben-Porat; M. Berwick; M. Robson; D. Kutler & Auerbach, A. D. "A Note on Competing Risks in Survival Data Analysis". In: British Journal of Cancer, 91 (2004), pp. 1229-1235.
}
\author{
Robin Denz, <robin.denz@rub.de>
}
