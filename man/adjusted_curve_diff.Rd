\name{adjusted_curve_diff}
\alias{adjusted_curve_diff}

\title{
Test if there is a difference between two Confounder-Adjusted Survival Curves or CIFs
}
\description{
This function implements a modified version of the Pepe and Flemming (1989) test for the difference of two adjusted survival curves or CIFs. In particular, the Null-Hypothesis is that the integral of the difference of the two curves in a specified time interval is equal to zero.
}
\usage{
adjusted_curve_diff(adjsurv, to, from=0, conf_level=0.95)
}

\arguments{
  \item{adjsurv}{
An \code{adjustedsurv} object created using the \code{adjustedsurv} function, or a \code{adjustedcif} object created using the \code{adjustedcif} function, with \code{bootstap=TRUE} in the original function call.
  }
    \item{to}{
A number specifying the right side of the time interval of interest. It has to be a value of time that can be read from both of the estimated survival curves or CIFs.
  }
  \item{from}{
A number specifying the left side of the time interval of interest. It has to be a value of time that can be read from both of the estimated survival curves or CIFs. It is set to 0 by default.
  }
  \item{conf_level}{
A number specifying the confidence level of the bootstrap confidence intervals.
  }
}
\details{

The \code{adjustedsurv} and \code{adjustedcif} functions with \code{bootstrap=TRUE} draw \code{n_boot} bootstrap samples and estimate the adjusted curves for each one. This function uses those estimates and calculates the integral of the difference between two curves (just the two when there are two curves and pairwise comparisons if there are > 2 curves) in the interval defined by \code{to} and \code{from}. If the curves are approximately equal, this quantity should be close to zero. The direct variance calculation of this is quantity is quite involved even in the non-adjusted case and has not been proposed for adjusted survival curves or adjusted CIFs yet. We can however use the distribution of the integrals over all bootstrap samples to approximate the variation. By shifting the bootstrap distribution to be centered around 0 we approximate the distribution of the integral under the Null-Hypothesis. The P-Value can then be calculated by taking the proportion of cases where the absolute of the integral observed in the actual curves is smaller or equal to the shifted bootstrap values.

Instead of relying on numerical integration, this function uses step-function interpolation, treating all supplied curves as step functions. The integral is then calculated using a simple sum of the squares defined by each step.

When the \code{adjustedsurv} or \code{adjustedcif} object was fitted using multiply imputed datasets, the tests are performed separately for each dataset. The estimates for the integral of the difference are combined using Rubins Rule. The confidence intervals for this quantity are calculated by pooling the bootstrap standard errors and recalculating the confidence interval using the normal approximation. The P-Values are also pooled using a method described in Licht (2010). It is recommended to check if the pooled P-Value is in agreement with the pooled confidence interval.

}
\value{
Returns a \code{curve_test} object. If there are exactly two treatments this list contains the following object:

\item{diff_curves}{A \code{data.frame} containing the difference curves used for calculating the integrals.
}
\item{diff_intergals}{A numeric vector containing the integrals of the difference for the estimated adjusted survival curves or CIFs.
}
\item{observed_diff_curve}{The curve of the difference between the two non-bootstrapped adjusted survival curves or CIFs.
}
\item{observed_diff_integral}{The integral of the curve in \code{observed_diff_curve}.
}
\item{p_value}{The P-Value for the modified Pepe-Fleming Test. See details.
}
\item{n_boot}{The number of bootstrap repetitions used.
}
\item{kind}{Whether survival curves or cumulative incidence functions where used.
}
\item{conf_int}{The percentile bootstrap confidence interval of the difference between the two curves.
}
\item{categorical}{Whether there are more than two treatments/groups.
}
\item{treat_labs}{The labels of all treatments/groups.
}
\item{method}{The adjustment method used in the original \code{adjustedsurv} or \code{adjustedcif} object.
}
\item{call}{The original function call.
}

If there are more than two treatment groups the object returned is a list of these objects with one list for each pairwise comparison.

If multiply imputed datasets where used, the object also includes a \code{mids_analyses} object, including a \code{curve_test} object for each imputed dataset. It also includes a \code{mids_p_values} object containing the separately estimated P-Values.

}
\references{
Margaret Sullivan Pepe and Thomas R. Fleming (1989). "Weighted Kaplan-Meier Statistics: A Class of Distance Tests for Censored Survival Data". In: Biometrics 45.2, pp. 497-507

Margaret Sullivan Pepe and Thomas R. Fleming (1991). "Weighted Kaplan-Meier Statistics: Large Sample and Optimality Considerations". In: Journal of the Royal Statistical Society: Series B 53.2, pp. 341-352

Nicholas I. Fisher and Peter Hall (1990). "On Bootstrap Hypothesis Testing". In: Australian Journal of Statistics 32.2, pp. 177-190

Florent Le Borgne, Bruno Giraudeau, Anne Héléne Querard, Magali Giral, and Yohann Foucher (2016). "Comparisons of the Performance of Different Statistical Tests for Time-To-Event Analysis with Confounding Factors: Practical Illustrations in Kidney Transplantation". In: Statistics in Medicine 35, pp. 1103-1116

Christine Licht (2010). "New Methods for Generating Significance Levels from Multiply-Imputed Data". PhD thesis. Otto-Friedrich-Universität Bamberg, Fakultät Sozial- und Wirtschaftswissenschaften
}
\author{
Robin Denz
}

\seealso{
\code{\link{plot.curve_test}}, \code{\link{adjustedsurv}}, \code{\link{adjustedcif}}
}
\examples{
\dontrun{
library(survival)
library(cmprsk)

#### Simple Survival Case with adjusted survival curves ####

# simulate some data as example
sim_dat <- sim_confounded_surv(n=500, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# estimate a cox-regression for the outcome
cox_mod <- coxph(Surv(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6 + group,
                 data=sim_dat, x=TRUE)


# use it to calculate adjusted survival curves with bootstrapping
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="direct",
                        outcome_model=cox_mod,
                        conf_int=FALSE,
                        bootstrap=TRUE,
                        n_boot=500)

# test the equality of both curves in the interval 0 to 1
adjtest <- adjusted_curve_diff(adjsurv, from=0, to=1)
print(adjtest)

#### Competing Risks Case with adjusted CIFs ####
library(riskRegression)
library(prodlim)

# simulate some data as example
sim_dat <- sim_confounded_crisk(n=500, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# estimate a cause-specific cox-regression for the outcome
csc_mod <- CSC(Hist(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6 + group,
                 data=sim_dat)


# use it to calculate adjusted CIFs for cause = 1 with bootstrapping
adjcif <- adjustedcif(data=sim_dat,
                      variable="group",
                      ev_time="time",
                      event="event",
                      method="direct",
                      outcome_model=csc_mod,
                      conf_int=FALSE,
                      bootstrap=TRUE,
                      n_boot=500,
                      cause=1)

# test the equality of both curves in the interval 0 to 1
adjtest <- adjusted_curve_diff(adjcif, from=0, to=1)
print(adjtest)

}
}
