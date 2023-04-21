
#' d_StoreyBH
#'
#' @param S_Y : score vector for the test set
#' @param S_X : score vector for the calibration set
#' @param alpha : significance level. Default level is set equal to 0.1
#' @param lambda : parameter involved in the computation of Storey estimator. Default value is set equal to 0.5
#'
#' @return An integer which is the \eqn{(1 − \alpha)}-confidence lower bound for the
#' number of true discoveries in Benjamini-Hochberg procedure with Storey estimator
#' for the proportion of true null hypotheses applied to conformal *p*-values.
#' The selection set, i.e. the set of hypothesis indices that we are
#' interested in is \eqn{[m]=:\{1,...,m\}} by default.
#' Then, Storey estimator is computed as
#' \deqn{\hat\pi_0 = \frac{1+\sum_{i=1}^m \mathbb{1}\{p_i>\lambda\}}{m(1-\lambda)}}
#' where \eqn{\lambda\in(0,1)} and \eqn{p_i} is the *p*-value related
#' to the null hypothesis \eqn{H_i, \hspace{2mm} i=1,\ldots,m}.

#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' d_StoreyBH(S_Y=Sy, S_X=Sx)
#'
#'
#'
d_StoreyBH = function(S_Y, S_X, alpha = 0.1, lambda=0.5){
  m = length(S_Y)
  n = length(S_X)
  pval = sort(sapply(1:m, function(i) (1+sum(S_X >= S_Y[i]))/(n+1)), decreasing=FALSE)
  pi0Sto = (1+sum(pval>lambda))/(m*(1-lambda))
  d =  sum(stats::p.adjust(pval,"BH")<=alpha/pi0Sto)
  return(d)
}





