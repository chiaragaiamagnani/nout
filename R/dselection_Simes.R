



#' dselection_Simes
#'
#' @param S_Y : score vector for the test set
#' @param S_X : score vector for the calibration set
#' @param S : vector of selected indices
#' @param alpha : significance level
#'
#' @return An integer which is the \eqn{(1 − \alpha)}-confidence lower bound for
#' the number of true discoveries in closed testing procedure using
#' Simes local test applied to conformal *p*-values. The selection set, i.e., the set of hypothesis
#' indices that we are interested can be chosen.
#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' dselection_Simes(S_Y=Sy, S_X=Sx, S = 5:15)
#'
dselection_Simes = function(S_Y, S_X, S, alpha = 0.1){
  S_Y.selection = S_Y[S]
  n = length(S_Y.selection)
  m = length(S_X)
  pval = sapply(1:n, function(i) (1+sum(S_X >= S_Y.selection[i]))/(m+1))
  d = hommel::discoveries(hommel::hommel(pval), alpha = alpha)

  # n = length(S_Y)
  # m = length(S_X)
  # pval = sapply(1:n, function(i) (1+sum(S_X >= S_Y[i]))/(m+1))
  # d = hommel::discoveries(hommel::hommel(pval), ix = selection, alpha = alpha)

  return(d)
}



