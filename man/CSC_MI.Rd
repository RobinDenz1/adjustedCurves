\name{CSC_MI}
\alias{CSC_MI}

\title{
Cause-Specific Cox Regression with Multiple Imputation
}
\description{
This function can be utilized to perform Cause-Specific Cox Regression on multiply imputed datasets.
}
\usage{
CSC_MI(mids, formula, ...)
}

\arguments{
  \item{mids}{
A \code{mids} object created using the \code{\link[mice]{mice}} function. This replaces the \code{data} argument in the original function call.
  }
  \item{formula}{
A formula object passed to the \code{\link[riskRegression]{CSC}} function in the \pkg{riskRegression} package.
  }
  \item{...}{
Other arguments which should be passed to the \code{\link[riskRegression]{CSC}} function in the \pkg{riskRegression} package.
  }
}
\details{

A small convenience function to perform CSC regression on multiply imputed data. It is simply a wrapper around the \code{\link[riskRegression]{CSC}} function from the \pkg{riskRegression} package, because the usual use of \code{with} is not supported directly. It returns a \code{mira} object, which can be passed to the \code{outcome_model} argument inside of the \code{\link[adjustedCurves]{adjustedcif}} function when needed. No \code{pool} method or other functionality is available.

}
\value{

A \code{mira} object containing the CSC regression for every imputed dataset.

}

\author{
Robin Denz
}

\seealso{
\code{\link{adjustedsurv}}, \code{\link[riskRegression]{CSC}}, \code{\link[mice]{mice}}
}
\examples{
# not run because it would be too slow
\donttest{
library(adjustedCurves)
library(survival)

if (requireNamespace("riskRegression") & requireNamespace("mice")) {
library(riskRegression)
library(mice)

# simulate some data as example
sim_dat <- sim_confounded_crisk(n=50, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# introduce random missingness in x1 as example
sim_dat$x1 <- ifelse(runif(n=50) < 0.5, sim_dat$x1, NA)

# perform multiple imputation
mids <- mice::mice(data=sim_dat, method="pmm", m=5, printFlag=0)

# use the function
csc_mods <- CSC_MI(mids=mids,
                   formula=Hist(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6 + group
                   )
}
}
}
