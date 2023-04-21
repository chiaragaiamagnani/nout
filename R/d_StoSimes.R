
#' d_StoreySimes
#'
#' @description It returns the lower bound for the number of true discoveries in closed testing procedure
#' using Simes local test with Storey estimator for the proportion of true
#' null hypotheses.
#'
#' @param S_Y : score vector of test observations
#' @param S_X : score vector of calibration observations
#' @param alpha : significance level of the local test. Default value is set equal to 0.1
#' @param lambda : parameter involved in the computation of Storey estimator. Default value is set equal to 0.5
#'
#'
#' @return An integer which is the \eqn{(1 − \alpha)}-confidence lower bound for the
#' number of true discoveries in closed testing procedure using Simes local
#' test with Storey estimator for the proportion of true null hypotheses applied to conformal *p*-values.
#' The selection set, i.e. the set of hypothesis indices that we are
#' interested in is \eqn{[m]=:\{1,...,m\}} by default.
#' Then, Storey estimator is computed as
#' \deqn{\hat\pi_0 = \frac{1+\sum_{i=1}^m \mathbb{1}\{p_i>\lambda\}}{m(1-\lambda)}}
#' where \eqn{\lambda\in(0,1)} and \eqn{p_i} is the *p*-value related
#' to the null hypothesis \eqn{H_i, \hspace{2mm} i=1,\ldots,m}.
#'
#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' d_StoreySimes(S_Y=Sy, S_X=Sx)
#'
#'
#'
d_StoreySimes = function(S_Y, S_X, alpha = 0.1, lambda=0.5){
  m = length(S_Y)
  n = length(S_X)
  pval = sort(sapply(1:m, function(i)
    (1+sum(S_X >= S_Y[i]))/(n+1)), decreasing=FALSE)

  simes.pval = sapply(1:m, function(i)
    min(pval[i:m]/seq(from=1, to=m-i+1, by=1)))

  # Building the levels of the Simes test with Storey estimator
  pi.not = sapply(1:m, function(i)
    (1+sum(pval[i:m]>lambda))/((m-i+1)*(1-lambda)))
  coeff = seq(from = m, to = 1, by = -1)
  thr = alpha/(coeff*pi.not)

  d = sum(cumsum(simes.pval <= thr) == 1:m)

  return(d)
}



