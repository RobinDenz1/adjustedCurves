\name{adjusted_curve_test}
\alias{adjusted_curve_test}

\title{
Test if there is a difference between two Confounder-Adjusted Survival Curves or CIFs
}
\description{
This function implements a modified version of the Pepe and Flemming (1989) test for the difference between two adjusted survival curves or CIFs. In particular, the Null-Hypothesis is that the integral of the difference of the two curves in a specified time interval is equal to zero.
}
\usage{
adjusted_curve_test(adj, to, from=0, conf_level=0.95,
                    interpolation="steps",
                    group_1=NULL, group_2=NULL)
}

\arguments{
  \item{adj}{
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
  \item{interpolation}{
Either \code{"steps"} (default) or \code{"linear"}. This parameter controls how interpolation is performed. If this argument is set to \code{"steps"}, the curves will be treated as step functions. If it is set to \code{"linear"}, the curves wil be treated as if there are straight lines between the point estimates instead. Points that lie between estimated points will be interpolated accordingly. Should usually be kept at \code{"steps"}. See Details.
  }
  \item{group_1}{
Optional argument to get one specific hypothesis test. This argument takes a single character string specifying one of the levels of the \code{variable} used in the original \code{adjustedsurv} or \code{adjustedcif} function call. This group will be subtracted from. For example if \code{group_1="A"} and \code{group_2="B"} the difference \code{A - B} will be used. If \code{NULL}, the order of the factor levels in the original \code{data} determines the test order. If not \code{NULL}, the \code{group_2} argument also needs to be specified. When these arguments are used, all other potential pairwise comparisons are ignored.
  }
  \item{group_2}{
Also a single character string specifying one of the levels of \code{variable}. This corresponds to the right side of the difference equation. See argument \code{group_2}.
  }
}
\details{

The \code{adjustedsurv} and \code{adjustedcif} functions with \code{bootstrap=TRUE} draw \code{n_boot} bootstrap samples and estimate the adjusted curves for each one. This function uses those estimates and calculates the integral of the difference between two curves in the interval defined by \code{to} and \code{from}. If the curves are approximately equal, this quantity should be close to zero. The direct variance calculation of this is quantity is quite involved even in the non-adjusted case and has not been proposed for adjusted survival curves or adjusted CIFs yet. We can however use the distribution of the integrals over all bootstrap samples to approximate the variation. By shifting the bootstrap distribution to be centered around 0 we approximate the distribution of the integral under the Null-Hypothesis. The p-value can then be calculated by taking the proportion of cases where the absolute of the integral observed in the actual curves is smaller or equal to the shifted bootstrap values.

The associated \code{print} and \code{summary} methods can be used to obtain a neat data.frame of the most important quantities. We also recommend checking the test assumptions using the \code{\link[=plot.curve_test]{plot}} method.

\strong{\emph{Pairwise Comparisons}}

When there are more than two survival curves or CIFs this function automatically performs pairwise comparisons between those. It is recommended to adjust the p-values obtained using this method for multiple testing. See \code{?p.adjust} for more information. If only one of the potential pairwise comparisons is of interest, the \code{group_1} and \code{group_2} arguments can be used to obtain only this specific one.

\strong{\emph{Multiple Imputation}}

When the \code{adjustedsurv} or \code{adjustedcif} object was fitted using multiply imputed datasets, the tests are performed separately for each dataset. The estimates for the integral of the difference are combined using Rubins Rule. The confidence intervals for this quantity are calculated by pooling the bootstrap standard errors and recalculating the confidence interval using the normal approximation. The p-values are also pooled using a method described in Licht (2010). It is recommended to check if the pooled p-value is in agreement with the pooled confidence interval.

\strong{\emph{Graphical Displays}}

To plot the curves of the differences directly, we recommend using the  \code{\link{plot_curve_diff}} function. Similar to the main plot functions, it has a lot of arguments to customize the plot. If the main goal is to check the assumptions, we recommend using the associated \code{\link[=plot.curve_test]{plot}} method instead.

\strong{\emph{Computational Details}}

Instead of relying on numerical integration, this function uses exact calculations. This is achieved by using either step-function interpolation (\code{interpolation="steps"}, the default) or linear interpolation (\code{interpolation="linear"}). In the former case, the integral is simply the sum of the area of the squares defined by the step function. In the second case, the integral is simply the sum of the area of the rectangles. Either way, there is no need for approximations. In some situations (for example when using parametric survival models with \code{method="direct"}), the curves are not step functions. In this case the \code{interpolation} argument should be set to \code{"linear"}.

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
\item{integral_se}{The bootstrap standard error of the difference integral.
}
\item{p_value}{The p-value for the modified Pepe-Fleming Test. See details.
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
\item{interpolation}{The interpolation method specified in the original \code{adjustedsurv} or \code{adjustedcif} object.
}
\item{call}{The original function call.
}

If there are more than two treatment groups the object returned is a list of these objects with one list for each pairwise comparison.

If multiply imputed datasets where used, the object also includes a \code{mids_analyses} object, including a \code{curve_test} object for each imputed dataset. It also includes a \code{mids_p_values} object containing the separately estimated p-values.

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
library(adjustedCurves)
library(survival)
library(cmprsk)

#### Simple Survival Case with adjusted survival curves ####

# simulate some data as example
set.seed(42)
sim_dat <- sim_confounded_surv(n=50, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# estimate a cox-regression for the outcome
cox_mod <- coxph(Surv(time, event) ~ x1 + x2 + x3 + x4 + x5 + x6 + group,
                 data=sim_dat, x=TRUE)

# use it to estimate adjusted survival curves with bootstrapping
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="direct",
                        outcome_model=cox_mod,
                        conf_int=FALSE,
                        bootstrap=TRUE,
                        n_boot=10) # n_boot should be much higher in reality

# test the equality of both curves in the interval 0 to 1
adjtest <- adjusted_curve_test(adjsurv, from=0, to=1)
print(adjtest)

#### Competing Risks Case with adjusted CIFs ####
library(riskRegression)
library(prodlim)

# simulate some data as example
sim_dat <- sim_confounded_crisk(n=50, max_t=1.2)
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
                      n_boot=10,
                      cause=1)

# test the equality of both curves in the interval 0 to 1
adjtest <- adjusted_curve_test(adjcif, from=0, to=1)
print(adjtest)
}
