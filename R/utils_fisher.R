# -------------------------------------------------------------------------------------- #
#  Implementation of Fisher's method (Bates et al. 2023) with asymptotical distribution  #
# -------------------------------------------------------------------------------------- #




#' stat.Fisher
#'
#' @description It computes the inverse of conformal *p*-values for each observations in the test set.
#' Then, they are transformed applying the logarithmic function.
#'
#' @param Z : pooled score vector with the first \eqn{m} components corresponding
#' to the calibration observations and the last \eqn{n} components corresponding
#' to the calibration observations
#' @param m : calibration sample size
#'
#'
#' @return Given the pooled score vector \eqn{Z=(X,Y)} where \eqn{X} is the calibration score
#' vector and \eqn{Y} is the test score vector, for each observation in the test sample
#' it returns the following quantity: \deqn{-log(p_i),}
#' where \eqn{p_i = \frac{1+\sum{j=1}^m\mathbb{1}\{X_j<Y_i\}}{m+1}} is the conformal
#' *p*-value of the \eqn{i}th observation in the test set.
#'
#'
stat.Fisher <- function(Z, m) {
  N = length(Z)
  n = N-m
  S_X = Z[1:m]
  S_Y = Z[(m+1):N]
  pval = sapply(1:n, function(i) (1+sum(S_X >= S_Y[i]))/(m+1))
  R = -2 * log(pval)
  return(R)
}


#' asymptotic.critical.Fisher
#'
#' @description It computes the \eqn{(1-\alpha)}-quantile of the adjusted Fisher test statistic
#' based on asymptotic Chi-squared approximation.

#' @param m : calibration sample size
#' @param n : test sample size
#' @param alpha : significance level. Default value is set equal to 0.1
#'
#'
#' @return It returns the \eqn{(1-\alpha)}-quantile of adjusted Fisher test statistic
#' based on asymptotic Chi-squared approximation
#' \deqn{\frac{-2\sum_{i=1}^n\log(p_i) + 2n\sqrt{\gamma+1}-1}{\sqrt{\gamma+1}} \overset{d}{\sim} \chi^2_{2n},}
#' where \eqn{\gamma\in(0,+\infty)} is chosen equal to \eqn{n/m}, \eqn{p_i} are the conformal *p*-values and
#' \eqn{\chi^2_{2n}} is the Chi-squared distribution with \eqn{2n} degrees of freedom.
#'
#'
asymptotic.critical.Fisher <- function(m, n, alpha=0.1) {
  gamma = n/m
  critical.value = sqrt(1+gamma) * stats::qchisq(p=1-alpha, df=2*n) - 2 * (sqrt(1+gamma)-1) * n
  return(critical.value)
}


#' asymptotic.pvalue.Fisher
#'
#' @description It computes the approximated *p*-value of the adjusted Fisher
#' test statistic based on the asymptotic Chi-squared approximation.
#'
#' @param m : calibration sample size
#' @param n : test sample size
#' @param T.obs : observed value of the test statistic
#'
#'
#' @return It returns the approximated *p*-value of the adjusted Fisher
#' test statistic based on the asymptotic Chi-squared approximation
#' \deqn{\frac{-2\sum_{i=1}^n\log(p_i) + 2n\sqrt{\gamma+1}-1}{\sqrt{\gamma+1}} \overset{d}{\sim} \chi^2_{2n},}
#' where \eqn{\gamma\in(0,+\infty)} is chosen equal to \eqn{n/m}, \eqn{p_i} are the conformal *p*-values and
#' \eqn{\chi^2_{2n}} is the Chi-squared distribution with \eqn{2n} degrees of freedom.
#'
#'
asymptotic.pvalue.Fisher <- function(m, n, T.obs) {
  gamma = n/m
  T.obs_shifted = (T.obs+2 * (sqrt(1+gamma)-1) * n)/ sqrt(1+gamma)
  p.value = stats::pchisq(q=T.obs_shifted, df=2*n, lower.tail = F)
  return(p.value)
}

