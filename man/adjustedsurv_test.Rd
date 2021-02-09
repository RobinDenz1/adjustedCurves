\name{adjustedsurv_test}
\alias{adjustedsurv_test}

\title{
Test for the Equality of Two Confounder-Adjusted Survival Curves
}
\description{
This function implements a modified version of the Pepe and Flemming (1989) test for the equality of two adjusted survival curves. In particular, the Null-Hypothesis is that the integral of the difference of the two curves in a specified time interval is not different from zero. The input of this function should be a \code{adjustedsurv} object created using the \code{adjustedsurv} function, with \code{bootstap=TRUE} in the original function call.
}
\usage{
adjustedsurv_test(adjsurv, from=0, to=Inf)
}

\arguments{
  \item{adjsurv}{
An \code{adjustedsurv} object created using the \code{adjustedsurv} function.
  }
  \item{from}{
A number specifying the left side of the time interval of interest. See details.
  }
  \item{to}{
A number specifying the right side of the time interval of interest. See details.
  }
}
\details{

To \code{adjustedsurv} function with \code{bootstrap=TRUE} draws \code{n_boot} bootstrap samples and estimates the adjusted survival curves for each one. This function uses those estimates and calculates the integral of the difference between two curves (just the two when there are two curves and pairwise comparisons if there are > 2 curves) in the interval defined by \code{to} and \code{from}. If the survival curves are approximately equal, this quantity should be close to zero. The direct variance calculation of this is quantity is quite involved even in the non-adjusted case (see ZITATE) and has not been proposed for adjusted survival curves yet. We can however use the distribution of the integrals over all bootstrap samples to approximate the variation. By shifting the bootstrap distribution to be centered around 0 we approximate the distribution of the integral under the Null-Hypothesis. The P-Value can then be calculated by taking the proportion of cases where the absolute of the integral observed in the actual curves is smaller or equal to the shifted bootstrap values.

This function has not been tested extensively yet and is still in the experimental stage. Use with caution.

}
\value{
Returns a \code{adjustedsurv_test} object containing the following object:

\item{diff_curves}{A \code{data.frame} containing the difference curves used for calculating the integrals.
}
\item{diff_intergals}{A numeric vector containing the integrals of the difference for the estimated adjusted survival curves.
}
\item{observed_diff_curve}{The curve of the difference between the two non-bootstrapped adjusted survival curves.
}
\item{obersved_diff_integral}{The integral of the curve in \code{observed_diff_curve}.
}
\item{p_value}{The P-Value for the modified Pepe-Fleming Test. See details.
}

}
\references{
Margaret Sullivan Pepe and Thomas R. Fleming. "Weighted Kaplan-Meier Statistics: A Class of Distance Tests for Censored Survival Data". In: Biometrics 45.2 (1989), pp. 497-507.

Margaret Sullivan Pepe and Thomas R. Fleming. "Weighted Kaplan-Meier Statistics: Large Sample and Optimality Considerations". In: Journal of the Royal Statistical Society: Series B 53.2 (1991), pp. 341–352.

Nicholas I. Fisher and Peter Hall. "On Bootstrap Hypothesis Testing". In: Australian Journal of Statistics 32.2 (1990), pp. 177-190.

Florent Le Borgne, Bruno Giraudeau, Anne Héléne Querard, Magali Giral, and Yohann Foucher. "Comparisons of the Performance of Different Statistical Tests for Time-To-Event Analysis with Confounding Factors: Practical Illustrations in Kidney Transplantation". In: Statistics in Medicine 35 (2016), pp. 1103-1116.
}
\author{
Robin Denz
}

\seealso{
\code{\link{print.adjustedsurv_test}}, \code{\link{adjustedsurv}}
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
