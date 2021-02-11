\name{surv_method_el}
\alias{surv_method_el}

\title{
Empirical Likelihood Estimation Survival Curves
}
\description{
This method utilizes the empirical likelihood estimation methodology to calculate adjusted survival curves. Neither outcome nor treatment assignment models are required for this function to work. Since the exact form of both mechanisms are left unspecified, it is more robust to model misspecification than IPTW or direct adjustment. See details for more information.
}
\usage{
surv_method_el(data, variable, ev_time, event, conf_int, conf_level=0.95,
               times, treatment_vars, moment="first", standardize=F, ...)
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
  \item{conf_int}{
A logical variable, indicating whether the asymptotic variances and confidence intervals of the survival probabilities should be calculated. Currently not available for this method. The argument is only included for consistency.
  }
  \item{conf_level}{
A number specifying the confidence level of the confidence intervals. Currently not used.
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
%%  ~~ If necessary, more details than the description above ~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
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
