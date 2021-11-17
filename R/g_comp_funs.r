#### For Adjusted CIFs ####

## fastcmprsk

# multiply betas with covariates and take sum
apply_betas <- function(x, betas) {
  return(sum(x * betas))
}

# using the same method as in predict.fccr, but vectorised to get the
# predictions in matrix form first
average_CIF_fccr <- function(object, newdata) {

  # apply betas
  eff <- apply(X=newdata, MARGIN=1, FUN=apply_betas, betas=object$coef)

  # same as in predict.fccr
  CIF.hat <- t(outer(X=exp(eff), Y=object$breslowJump[, 2], FUN="*"))
  CIF.hat <- apply(X=CIF.hat, MARGIN=2, FUN=cumsum)
  CIF.hat <- apply(X=CIF.hat, MARGIN=2, FUN=function(x){1 - exp(-x)})

  # take average at each t
  CIF_z <- apply(X=CIF.hat, MARGIN=1, FUN=mean)

  return(CIF_z)
}

##
