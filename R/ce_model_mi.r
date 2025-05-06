
ce_model_mi <- function(mids, formula, cause=1, model, ...) {

  args <- list(...)
  if (!is.null(args$data)) {
    stop("The 'data' argument cannot be used here. It is replaced by the",
         "'mids' argument.")
  }

  imp_long <- mice::complete(mids, action="long", include=FALSE)
  outc_mod <- list()
  for (i in seq_len(max(imp_long$.imp))) {
    if (model=="CSC") {
      mod <- riskRegression::CSC(formula,
                                 data=imp_long[imp_long$.imp==i, ],
                                 ...)
    } else if (model=="FGR") {
      mod <- riskRegression::FGR(formula,
                                 data=imp_long[imp_long$.imp==i, ],
                                 cause=cause,
                                 ...)
    }
    mod$call$formula <- formula
    outc_mod[[i]] <- mod
  }

  mira_obj <- list(call=match.call(),
                   call1=mids$call,
                   n_mis=mids$nmis,
                   analyses=outc_mod)
  class(mira_obj) <- c("mira", "matrix")
  return(mira_obj)

}

## fit a Fine & Gray regression with multiply imputed data
#' @export
FGR_MI <- function(mids, formula, cause=1, ...) {
  return(ce_model_mi(mids=mids, formula=formula, cause=cause,
                     model="FGR", ...))
}

## fit a cause-specific cox regression with multiply imputed data
#' @export
CSC_MI <- function(mids, formula, ...) {
  return(ce_model_mi(mids=mids, formula=formula, model="CSC", ...))
}
