



#' d.S_simes
#'
#' @description It returns the lower bound for the number of true discoveries in closed testing procedure
#' using Simes local test for any selection in the test index set.
#'
#'
#' @param S_Y : score vector for the test set
#' @param S_X : score vector for the calibration set
#' @param S : vector of selected indices in the test set
#' @param alpha : significance level. Default values is set equal to 0.1
#'
#' @return An integer which is the \eqn{(1 − \alpha)}-confidence lower bound for
#' the number of true discoveries in closed testing procedure using
#' Simes local test applied to conformal *p*-values.
#' The selection set can be any arbitrary subset in the index set of test observations
#' according to which test units we want to test.
#'
#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' d.S_simes(S_Y=Sy, S_X=Sx, S = 5:15)
#'
d.S_simes = function(S_Y, S_X, S, alpha = 0.1){

  n = length(S_Y)
  m = length(S_X)
  pval = sapply(1:n, function(i) (1+sum(S_X >= S_Y[i]))/(m+1))
  d = hommel::discoveries(hommel::hommel(pval), ix = S, alpha = alpha)

  return(d)
}



