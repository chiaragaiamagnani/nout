
#' d_benjhoch
#'
#' @param S_Y : score vector for the test set
#' @param S_X : score vector for the calibration set
#' @param alpha : significance level. Default level is set equal to 0.1
#'
#' @return An integer which is the \eqn{(1 − \alpha)}-confidence lower bound for
#' the number of true discoveries in Benjamini-Hochberg procedure applied to.
#' The selection set, i.e. the set of hypothesis
#' indices that we are interested in is \eqn{[m]=:\{1,...,m\}} by default.
#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' d_benjhoch(S_Y=Sy, S_X=Sx)
#'
#'
d_benjhoch = function(S_Y, S_X, alpha = 0.1){
  m = length(S_Y)
  n = length(S_X)
  pval = sapply(1:m, function(i) (1+sum(S_X >= S_Y[i]))/(n+1))
  d =  sum(stats::p.adjust(pval,"BH")<=alpha)
  return(d)
}






