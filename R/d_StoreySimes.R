
#' d_storeysimes
#'
#' @description It returns the lower bound for the number of true discoveries in closed testing procedure
#' using Simes local test with Storey estimator for the proportion of true null hypotheses.
#'
#' @param S_Y : score vector of test observations
#' @param S_X : score vector of calibration observations
#' @param alpha : significance level of the local test. Default value is set equal to 0.1
#' @param lambda : parameter involved in the computation of Storey estimator. Default value is set equal to 0.5
#'
#'
#' @return An integer which is the \eqn{(1 − \alpha)}-confidence lower bound for the
#' number of true discoveries in closed testing procedure using Simes local
#' test with Storey's estimator for the proportion of true null hypotheses applied to conformal *p*-values.
#' The selection set is trivial, i.e., we are interested in testing all the observations in the test set by default.
#' Then, Storey estimator is computed as
#' \deqn{\hat\pi_0 = \frac{1+\sum_{i=1}^n \mathbb{1}\{p_i>\lambda\}}{n(1-\lambda)}}
#' where \eqn{n} is the test sample size, \eqn{\lambda\in(0,1)} is a tuning parameter
#' and \eqn{p_i} is the *p*-value corresponding to the \eqn{i}th hypothesis in the test set.
#'
#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' d_storeysimes(S_Y=Sy, S_X=Sx)
#'
#'
#'
d_storeysimes = function(S_Y, S_X, alpha = 0.1, lambda=0.5){
  m = length(S_Y)
  n = length(S_X)
  pval = sort(sapply(1:m, function(i)
    (1+sum(S_X >= S_Y[i]))/(n+1)), decreasing=FALSE)

  simes.pval = sapply(1:m, function(i)
    min(pval[i:m]/seq(from=1, to=m-i+1, by=1)))

  # Building the levels of Simes test with Storey estimator
  # pi.not = sapply(1:m, function(i)
  #   (1+sum(pval[i:m]>lambda))/((m-i+1)*(1-lambda)))

  # Building the levels of Simes test with Storey estimator.
  # Storey estimator will be used in the closed testing procedure
  # in every levels except for lowest ones, when the set of
  # considered pvalues has cardinality less than or equal to 2.

  pi.not.highlevels = sapply(1:(m-2), function(i)
      (1+sum(pval[i:m]>lambda))/((m-i+1)*(1-lambda)))
  pi.not = c(pi.not.highlevels,1,1)
  coeff = seq(from = m, to = 1, by = -1)
  thr = alpha/(coeff*pi.not)

  d = sum(cumsum(simes.pval <= thr) == 1:m)

  return(list("d"=d, "pi.not"=pi.not[1]))
}



