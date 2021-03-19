\name{surv_method_el}
\alias{surv_method_el}

\title{
Empirical Likelihood Estimation Survival Curves
}
\description{
This method utilizes the empirical likelihood estimation methodology to calculate adjusted survival curves. Neither outcome nor treatment assignment models are required for this function to work. Since the exact form of both mechanisms are left unspecified, it is more robust to model misspecification than IPTW or direct adjustment. See details for more information.
}
\usage{
surv_method_el(data, variable, ev_time, event, times, treatment_vars,
               moment="first", standardize=F, ...)
}

\arguments{
  \item{data}{
A \code{data.frame} object containing the needed time-to-event data in standard format.
}
  \item{variable}{
A character string specifying the variable by which the survival curves should be grouped. Must be a valid columname of \code{data}.
}
  \item{ev_time}{
A character string specifying the variable indicating the time-to-event or time-to-censoring. Must be a valid columname of \code{data}.
}
  \item{event}{
A character string specifying the binary event indicator. Must be a valid columname of \code{data}.
}
  \item{times}{
Method used in \code{WeightIt} function call. Ignored if \code{treatment_model} is not a formula object.
}
  \item{treatment_vars}{
A character vector of column names specifying variables to be used as covariates in the empirical likelihood estimation. See details and examples.
}
  \item{moment}{
A character string specifying which moment to adjust for. See details.
}
  \item{standardize}{
A logical variable indicating whether the \code{treatment_vars} variables should be standardized. See details.
}
  \item{...}{
Further arguments passed to the \code{el.est} function inside of the \code{adjKMtest} package.
}
}
\details{

\itemize{
\item{\strong{Type of Adjustment:} Requires a character vector of variable names used to balance the distribution of covariates (treatment assignment mechanism)}
\item{\strong{Doubly-Robust:} Estimates are not Doubly-Robust.}
\item{\strong{Categorical groups:} Only binary treatments are allowed. Column specified by \code{variable} must be an integer with only 0 and 1 as values.}
\item{\strong{Approximate Variance:} Calculations to approximate the variance and confidence intervals are not available. The user has to use bootstrapping if those are desired.}
\item{\strong{Allowed Time Values:} Allows both continuous and integer time.}
\item{\strong{Bounded Estimates:} Estimates are guaranteed to be bounded in the 0 to 1 probability range.}
\item{\strong{Monotone Function:} Estimates are guaranteed to be monotone.}
\item{\strong{Dependencies:} This method relies on the \code{adjKMtest} package. It can be installed using the \code{install_github} function in the \code{devtools} R-Package (See \code{?install_github}) from the official repository: <https://github.com/kimihua1995/adjKMtest>}
}

}
\value{
Returns a \code{data.frame} object containing the estimated adjusted survival probabilities for some points in time for each level of \code{variable}.
}
\references{
Xiaofei Wang, Fangfang Bai, Herbert Pang, and Stephen L. George. "Bias-Adjusted Kaplan-Meier Survival Curves for Marginal Treatment Effect in Observational Studies". In: Journal of Biopharmaceutical Statistics 29.4 (2019), pp. 592-605.

Art B. Owen. Empirical Likelihood. Boca Raton: CRC Press, 2001.
}
\author{
The wrapper function was written by Robin Denz, the \code{el.est} function (which this wrapper is build around) was written by Kaiyuan "Kimi" Hua. See \code{?el.est} for more details.
}

\seealso{
\code{\link[adjKMtest]{el.est}}
}
\examples{
\dontrun{
# simulate some data as example (needs a binary integer "variable")
sim_dat <- sim_confounded_surv(n=500, max_t=1.2)

# calculate adjusted survival curves
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="el",
                        treatment_vars=c("x1", "x2", "x3", "x4", "x5", "x6"),
                        moment="first")

# plot the curves
plot(adjsurv)
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
