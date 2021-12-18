# Copyright (C) 2021  Robin Denz
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

## perform a cause-specific cox regression with multiply imputed data
#' @export
FGR_MI <- function(mids, formula, cause=1, ...) {

  args <- list(...)
  if (!is.null(args$data)) {
    stop("The 'data' argument cannot be used here. It is replaced by the",
         "'mids' argument.")
  }

  imp_long <- mice::complete(mids, action="long", include=FALSE)
  outc_mod <- list()
  for (i in seq_len(max(imp_long$.imp))) {
    mod <- riskRegression::FGR(formula,
                               data=imp_long[imp_long$.imp==i, ],
                               cause=cause,
                               ...)
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
