\name{adjustedCurves-package}
\alias{adjustedCurves-package}

\title{
Confounder-Adjusted Survival Curves and Cumulative Incidence Functions
}
\description{
This package includes various methods to calculate confounder-adjusted survival curves (single event) and confounder adjusted cumulative incidence functions (possibly with multiple competing events). It provides \code{plot} functions to easily produce highly customizable and publication-ready graphics. It also allows the user to easily calculate relevant statistics, such as confidence intervals, p-values, and adjusted restricted mean survival time estimates. Multiple Imputation is directly supported. The aim is to unite the sometimes very different adjustments methods under one consistent framework and to make the usage of these methods and the calculation of associated statistics as easy as possible for the user, while still providing substantial functionality.

The design of this package is based on the design of the \pkg{WeightIt} package. It includes two main functions: \code{\link[adjustedCurves]{adjustedsurv}} and \code{\link[adjustedCurves]{adjustedcif}}. Every implemented adjustment method has their own documentation page including a small description, code examples, and relevant literature references. The typical workflow using this package is as follows (1) calculate confounder-adjusted curves (survival curves or CIFs) using either \code{adjustedsurv} or \code{adjustedcif}, (2) plot those using the S3 \code{plot} method and sometimes (3) calculate further statistics using \code{adjusted_rmst} or \code{test_curve_equality}.

Bug reports, suggestions and feature requests are highly welcome. Please file an issue on the official github page or contact the author directly using the supplied e-mail address.
}
\author{
Robin Denz, <robin.denz@rub.de>
}
