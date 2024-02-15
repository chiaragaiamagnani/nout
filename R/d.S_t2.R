



#' d.S_t2
#'
#'
#' @description For any aribitrary selection subset of the index test set, the function returns a confidence
#' lower bound for the number of true discoveries provided by closed testing procedure using Wilcoxon-Mann-Whitney local test.
#'
#' @param S_X : score vector for the calibration set
#' @param S_Y : score vector for the test set
#' @param S : vector of selected indices in the test set
#' @param alpha : significance level. Default level is set equal to 0.1
#' @param n.exact : if \eqn{min\{m,n\}\leq n.exact} the critical values of the Wilcoxon-Mann-Whitney
#' statistic are exactly computed using ***qwilcox*** function. Default value is set equal to 10
#'
#'
#' @return An integer which is the \eqn{(1 − \alpha)}-confidence lower bound for
#' the number of true discoveries in closed testing procedure using
#' Wilcoxon-Mann-Whitney local test that rejects the global null hypothesis at level \eqn{\alpha}
#' when the test statistic \deqn{\sum_{i=1}^n\sum_{j=1}^m \mathbb{1}\{X_j<Y_i\}} is greater than
#' the critical value corresponding to the significance level \eqn{\alpha}.
#' The selection set can be any arbitrary subset in the index test set
#' according to which units we want to test.
#'
#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' d.S_t2(S_Y=Sy, S_X=Sx, S=c(1:3, 21), alpha=0.1)
#'
#'
d.S_t2 = function(S_X, S_Y, S, alpha=0.1, n.exact=10){

  n = length(S_Y)
  m = length(S_X)

  # statistics = ranks of S_Y[i] in (S_X,S_Y[i])
  g = sapply(1:n, function(i) sum(S_Y[i]>S_X))

  # Compute critical values
  if(min(n,m)<=n.exact){
    if(n.exact>20){
      stop("Execution stopped: exact computation of critical values is computationally too demanding.")
    }
    cvs = nout::exact.crit.t2(m=m, n=n, alpha=alpha)$crit.val
  }

  if(min(n,m)>n.exact){
    cvs = nout::exact.crit.t2(m=m, n=n, alpha=alpha)$crit.val
  }

  res = sumSome::sumStatsPar(g = g, S = S, alpha = alpha, cvs = cvs)

  return(list("d.S_t2"=res$TD, "S"=S))

}


