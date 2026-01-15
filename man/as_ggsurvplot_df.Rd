\name{as_ggsurvplot_df}
\alias{as_ggsurvplot_df}

\title{
Extract a \code{data.frame} containing the estimated survival curves from a \code{adjustedsurv} object
}
\description{
A small convenience function to extract the most important quantities from an \code{adjustedsurv} object. The resulting \code{data.frame} is structured according to the format required by the \code{ggsurvplot_df} function of the \pkg{survminer} package, making it easy to use the \code{ggsurvplot_df} function.
}
\usage{
as_ggsurvplot_df(adjsurv)
}
\arguments{
  \item{adjsurv}{
An object of class \code{adjustedsurv} created by the \code{\link{adjustedsurv}} function.
  }
}
\value{
Returns a \code{data.frame} containing the required information, extracted from the \code{adjustedsurv} object.
}
\author{
Robin Denz
}
\details{
If the \code{adjsurv} object was created using multiply imputed data, all values correspond to the pooled estimates, including numbers at risk and numbers of events. Further, if a weighted analysis was performed, the weights are used to derive the numbers at risk and numbers of events.
}
\seealso{
\code{\link{adjustedsurv}}, \code{\link{plot.adjustedsurv}}
}
\examples{
library(adjustedCurves)
library(survival)

set.seed(42)

# simulate some example data
sim_dat <- sim_confounded_surv(n=50, max_t=1.2)
sim_dat$group <- as.factor(sim_dat$group)

# treatment assignment model
glm_mod <- glm(group ~ x2 + x3 + x5 + x6, data=sim_dat, family="binomial")

# estimate some adjusted survival curves
adjsurv <- adjustedsurv(data=sim_dat,
                        variable="group",
                        ev_time="time",
                        event="event",
                        method="iptw_km",
                        treatment_model=glm_mod,
                        conf_int=TRUE,
                        bootstrap=FALSE)

# extract info
df <- as_ggsurvplot_df(adjsurv)

# not run here to avoid dependency on survminer
if (interactive()) {
# plot using survminer, requires the 'survminer' package
ggsurvplot_df(df)
}
}
