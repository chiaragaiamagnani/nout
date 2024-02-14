
#' discoveries_StoreyBH
#'
#' @description Applies Benjamini-Hochberg procedure to conformal *p*-values using
#' Storey estimator for the proportion of true null hypotheses.
#'
#' @param S_Y : score vector for the test set
#' @param S_X : score vector for the calibration set
#' @param alpha : significance level. Default level is set equal to 0.1
#' @param lambda : parameter involved in the computation of Storey estimator. Default value is set equal to 0.5
#'
#' @return A vector corresponding to indices of test observations which are rejected
#' by Benjamini-Hochberg procedure applied to conformal *p*-values and using Storey
#' estimator for the proportion of true null hypoteses.
#' The selection set is trivial, i.e., we are interested in testing all the observations in the test set by default.
#' Storey estimator is computed as
#' \deqn{\hat\pi_0 = \frac{1+\sum_{i=1}^n \mathbb{1}\{p_i>\lambda\}}{m(1-\lambda)}}
#' where \eqn{\lambda\in(0,1)} is a tuning parameter, \eqn{m} and \eqn{n} are the size of the calibration set and of the test set
#' respectively and \eqn{p_i} is the *p*-value related to the *i*th observation in the test set.

#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' discoveries_StoreyBH(S_Y=Sy, S_X=Sx)
#'
#'
#'
discoveries_StoreyBH = function(S_Y, S_X, alpha = 0.1, lambda=0.5){
  m = length(S_Y)
  n = length(S_X)
  pval = sort(sapply(1:m, function(i) (1+sum(S_X >= S_Y[i]))/(n+1)), decreasing=FALSE)
  pi0Sto = (1+sum(pval>lambda))/(m*(1-lambda))
  index.lesslambda = which(pval<lambda)
  D =  which(stats::p.adjust(pval,"BH")[index.lesslambda]<=alpha/pi0Sto)
  return(D)
}





