
## calculates the confidence interval of the ratio of two survival
## probabilities given the probabilities and their standard errors
# NOTE: This is a simplified version of the usual equation because we
#       can assume that a and b are independent (cov(a, b) = 0)
fieller_ratio_ci <- function(a, b, a_se, b_se, conf_level=0.95) {

  ratio <- a/b
  a_var_se <- a_se^2
  b_var_se <- b_se^2

  alpha <- 1 - conf_level
  z <- stats::qnorm(1 - alpha/2)
  g <- (z^2) * b_var_se/b^2
  C <- sqrt(a_var_se + ratio^2 * b_var_se - g*a_var_se)

  ci_lower <- (1/(1-g)) * (ratio - z/b * C)
  ci_upper <- (1/(1-g)) * (ratio + z/b * C)

  out <- list(ratio=ratio, ci_lower=ci_lower, ci_upper=ci_upper)

  return(out)
}

## calculates a p-value for a one-sample ratio test
fieller_p_val <- function(a, b, a_se, b_se, ratio=1) {
  a_var_se <- a_se^2
  b_var_se <- b_se^2

  t_stat <- (a - ratio * b) / sqrt(a_var_se + ratio^2 * b_var_se)
  p_val <- 2 * stats::pnorm(abs(t_stat), lower.tail=FALSE)

  return(p_val)
}
