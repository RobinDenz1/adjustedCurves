
## calculates the confidence interval of the ratio of two survival
## probabilities given the probabilities, their standard errors and
## their covariance
# NOTE: With ab_cov=0 (default) this is the simplified version of the
#       usual equation, which assumes that a and b are independent
fieller_ratio_ci <- function(a, b, a_se, b_se, conf_level=0.95, ab_cov=0) {

  ratio <- a/b
  a_var_se <- a_se^2
  b_var_se <- b_se^2

  alpha <- 1 - conf_level
  z <- stats::qnorm(1 - alpha/2)
  g <- (z^2) * b_var_se/b^2
  # cov(a, b) / var(b), set to 0 if b has no variance
  k <- ifelse(b_var_se==0, 0, ab_cov/b_var_se)
  C <- sqrt(a_var_se - 2*ratio*ab_cov + ratio^2 * b_var_se -
              g*(a_var_se - k*ab_cov))

  ci_lower <- (1/(1-g)) * (ratio - g*k - z/b * C)
  ci_upper <- (1/(1-g)) * (ratio - g*k + z/b * C)

  out <- list(ratio=ratio, ci_lower=ci_lower, ci_upper=ci_upper)

  return(out)
}

## calculates a p-value for a one-sample ratio test
fieller_p_val <- function(a, b, a_se, b_se, ratio=1, ab_cov=0) {
  a_var_se <- a_se^2
  b_var_se <- b_se^2

  t_stat <- (a - ratio * b) / sqrt(a_var_se + ratio^2 * b_var_se -
                                     2 * ratio * ab_cov)
  p_val <- 2 * stats::pnorm(abs(t_stat), lower.tail=FALSE)

  return(p_val)
}
