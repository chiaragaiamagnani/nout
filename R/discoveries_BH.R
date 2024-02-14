
#' discoveries_BH
#'
#' @description Applies Benjamini-Hochberg procedure to conformal *p*-values.
#'
#' @param S_Y : score vector for the test set
#' @param S_X : score vector for the calibration set
#' @param alpha : significance level. Default level is set equal to 0.1
#'
#' @return A vector corresponding to indices of test observations which are rejected
#' by Benjamini-Hochberg procedure applied to conformal *p*-values.
#' The selection set is trivial, i.e., we are interested in testing all the observations in the test set by default.
#'
#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' discoveries_BH(S_Y=Sy, S_X=Sx)
#'
#'
discoveries_BH = function(S_Y, S_X, alpha = 0.1){
  m = length(S_Y)
  n = length(S_X)
  pval = sapply(1:m, function(i) (1+sum(S_X >= S_Y[i]))/(n+1))
  D =  which(stats::p.adjust(pval,"BH")<=alpha)
  return(D)
}






