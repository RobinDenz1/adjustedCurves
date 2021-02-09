\name{sim_confounded_surv}
\alias{sim_confounded_surv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate Survival Data with Confounders
}
\description{
A function to simulate time-to-event data with one or multiple confounders. The user can specify both the relationship between the covariates and the survival time and the relationship between the covariates and the treatment assignment probability. Random censoring based on a custom function may also be introduced. Can be used for simulation studies or to showcase the useage of the adjusted survival curve methodology presented in this package.
}
\usage{
sim_confounded_surv(n=500, lcovars=NULL, outcome_betas=NULL,
                    group_beta=-1, surv_dist="weibull",
                    gamma=1.8, lambda=2, treatment_betas=NULL,
                    intercept=-0.5, gtol=0.001,
                    cens_fun=function(n){stats::rweibull(n, 1, 2)},
                    cens_args=list(), max_t=Inf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
  An integer specifying the sample size of the simulated data set.
  }
  \item{lcovars}{
  A named list to specify covariates. Each list element should be a vector containing information on the desired covariate distribution. See details.
  }
  \item{outcome_betas}{
  A named numeric vector of beta coefficients for the time-to-event outcome.
  }
  \item{group_beta}{
  A number specifying the beta coefficient of the grouping variable on the survival time.
  }
  \item{surv_dist}{
  A character string denoting the distribution used in the simulation of the survival time. See details.
  }
  \item{gamma}{
  A numeric parameter for the simulation of the survival time. See details.
  }
  \item{lambda}{
  A numeric parameter for the simulation of the survival time. See details.
  }
  \item{treatment_betas}{
  A named numeric vector of beta coefficients for the treatment assignment model.
  }
  \item{intercept}{
  The intercept of the treatment assignment model.
  }
  \item{gtol}{
  Tolerance at which estimated treatment assignment probabilities are truncated.
  }
  \item{cens_fun}{
  A function to generate censoring times.
  }
  \item{cens_args}{
  Arguments passed to \code{cens_fun}.
  }
  \item{max_t}{
  A number denoting the maximum follow-up time. Every event time bigger than this threshhold are censored.
  }
}
\details{

}
\value{
Returns a \code{data.frame} object containing the simulated covariates, the event indicator ("cause"), the survival/censoring time ("time") and the group variable ("group").
}

\references{
Hier mein eigenes Paper

Ralf Bender, Thomas Augustin, and Maria Blettner. "Generating Survival Times to Simulate Cox Proportional Hazards Models". In: Statistics in Medicine 24.11 (2005), pp. 1713-1723.

Arthur Chatton, Florent Le Borgne, Clemence Leyrat, and Yohann Foucher. G-Computation and Inverse Probability Weighting for Time-To-Event Outcomes: A Comparative Study. arXiv:2006.16859v1. (2020).
}
\author{
Robin Denz
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
