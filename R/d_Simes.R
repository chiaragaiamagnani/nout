#' d_Simes
#'
#' @description It returns the lower bound for the number of true discoveries in closed testing procedure
#' using Simes local test.
#'
#' @param S_Y : score vector of test observations
#' @param S_X : score vector of calibration observations
#' @param alpha : significance level of the local test. Default value is set equal to 0.1.
#'
#' @return An integer which is the \eqn{(1 − \alpha)}-confidence lower bound for the
#' number of true discoveries in closed testing procedure using Simes local
#' test applied to conformal *p*-values.
#' The selection set is trivial, i.e., we are interested in testing all the observations in the test set by default.
#'
#'
#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' d_Simes(S_Y=Sy, S_X=Sx)
#'
#'
d_Simes = function(S_Y, S_X, alpha = 0.1){
  n = length(S_Y)
  m = length(S_X)
  pval = sapply(1:n, function(i) (1+sum(S_X >= S_Y[i]))/(m+1))
  d = hommel::discoveries(hommel::hommel(pval), alpha = alpha)
  return(d)
}


