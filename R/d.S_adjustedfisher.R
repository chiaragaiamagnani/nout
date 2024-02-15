



#' d.S_adjustedfisher
#'
#' @description For any aribitrary selection subset of the index test set, the function returns a confidence
#' lower bound for the number of true discoveries provided by closed testing procedure using the adjusted Fisher
#' local test applied to conformal *p*-values.
#'
#' @param S_X : score vector for the calibration set
#' @param S_Y : score vector for the test set
#' @param S : vector of selected indices in the test set
#' @param alpha : significance level. Default level is set equal to 0.1
#' @param n.exact : if \eqn{min\{m,n\}\leq n.exact} the critical values of the adjusted Fisher test
#' statistic are computed via permutation. Default value is set equal to 10
#' @param B : number of iterations to approximate the adjusted Fisher statistic distribution via permutation
#' @param seed : seed to ensure reproducible results
#'
#' @return An integer which is the \eqn{(1 − \alpha)}-confidence lower bound for
#' the number of true discoveries in closed testing procedure using
#' the adjusted Fisher local test applied to conformal *p*-values. The test rejects
#' the global null hypothesis at level \eqn{\alpha} when
#' \deqn{\frac{-2\sum_{i=1}^n\log(p_i)+2n(\sqrt{1+\gamma}-1)}{\sqrt{1+\gamma}\geq \chi^2_{2n}(1-\alpha)},}
#' where *n* is the size of the test set and \eqn{\chi^2_{2n}(1-\alpha)} is the \eqn{(1-\alpha)} quantile of the chi-squared
#' distribution with \eqn{2n} degree of freedom.
#' The selection set can be any arbitrary subset in the index test set
#' according to which units we want to test.
#'
#' @references
#' Bates S., Candès E., Lei L., Romano Y. and Sesia M. (2023). Testing for outliers with conformal *p*-values. Annals of Statistics, doi: 10.1214/22-AOS2244
#'
#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' d.S_adjustedfisher(S_Y=Sy, S_X=Sx, S=c(1:3, 21), alpha=0.1)
#'
#'
#'
d.S_adjustedfisher = function(S_X, S_Y, S, alpha=0.1, n.exact=10, B=10^3, seed=123){
  n = length(S_Y)
  m = length(S_X)

  pval = sapply(1:n, function(i) (1+sum(S_X >= S_Y[i]))/(m+1))
  g = -2 * log(pval) # statistics

  # Compute critical values
  if(min(n,m)>n.exact){

    # Vector of critical values: the first entry is the critical value referred to
    # the first (up) level in the closed testing (where we test the global null) and
    # the last entry is the critical value referred to the last (bottom) level
    # in the closed testing (where we test the elementary hypotheses).
    # crit = (chi2_{2m}(1-alpha),chi2_{2(m-1)}(1-alpha),..., chi2_{2}(1-alpha))
    cvs = stats::qchisq(p = alpha, df = 2 * seq(n), lower.tail = FALSE)
  }

  if(min(n,m)<=n.exact){
    cvs = sapply(1:n, function(h) nout::perm.crit.adjustedfisher(S_X=S_X, S_Y=S_Y[1:h], B=B, alpha=alpha, seed=seed)$crit.val)
  }

  res = sumSome::sumStatsPar(g = g, S = S, alpha = alpha, cvs = cvs)

  return(list("d.S_adjustedfisher"=res$TD, "S"=S))

}



