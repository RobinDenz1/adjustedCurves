\name{CSC_MI}
\alias{CSC_MI}

\title{
Cause-Specific Cox Regression with multiple imputation
}
\description{
This function can be utilized to perform Cause-Specific Cox Regression on multiply imputed datasets.
}
\usage{
CSC_MI(mids, formula, ...)
}

\arguments{
  \item{mids}{
A \code{mids} object created using the \code{mice} function. This replaces the \code{data} argument in the original function call.
  }
  \item{formula}{
A formula object passed to the \code{CSC} function in the \code{riskRegression} package.
  }
  \item{...}{
Other arguments which should be passed to the \code{CSC} function in the \code{riskRegression} package.
  }
}
\details{

A small convenience function to perform CSC regression on multiply imputed data. It is simply a wrapper around the \code{CSC} function from the \code{riskRegression} package, because the usual use of \code{with} is not supported directly. It returns a \code{mira} object, which can be passed to the \code{outcome_model} argument inside of the \code{adjustedcif} function when needed. No \code{pool} method or other functionality is available.

}
\value{

A \code{mira} object containing the CSC regression for every imputed dataset.

}

\author{
Robin Denz
}

\seealso{
\code{\link{adjustedsurv}}
}
\examples{
\dontrun{
library(survival)
library(riskRegression)
library(mice)

# simulate some data as example
sim_dat <- sim_confounded_crisk(n=500, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# introduce random missingness in x1 as example
sim_dat$x1 <- ifelse(runif(n=500) < 0.5, sim_dat$x1, NA)

# perform multiple imputation
mids <- mice::mice(data=sim_dat, method="pmm", m=5)

# use the function
csc_mods <- CSC_MI(mids=mids,
                   formula=Hist(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6 + group
                   )


}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
