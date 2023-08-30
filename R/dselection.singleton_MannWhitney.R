
#' dselection.singleton_MannWhitney
#'
#' @description It returns the lower bound for the number of true discoveries in closed testing procedure
#' using Wilcoxon-Mann-Whitney local test.
#'
#' @param S_Y : score vector for the test set
#' @param S_X : score vector for the calibration set
#' @param S : selected set of indices. Must be a singleton.
#' @param alpha : significance level
#'
#' @return An integer which is the \eqn{(1 − \alpha)}-confidence lower bound for
#' the number of true discoveries in closed testing procedure using
#' Wilcoxon-Mann-Whitney local test applied to conformal *p*-values.
#' The selection set S, i.e. the set of hypothesis
#' indices that we are interested in, must be a singleton.
#'
#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' crit = nout::critWMW(m=length(Sx), n=length(Sy), alpha=0.1)
#' dselection.singleton_MannWhitney(S_Y=Sy, S_X=Sx, S=3, alpha=0.1)
#'
dselection.singleton_MannWhitney = function(S_Y,S_X,S,alpha){

  n = length(S_Y)
  m = length(S_X)
  s = length(S)

  if(s>1)
    stop("Error: S must be e singleton.")
  if(s==0)
    stop("Error: S is empty.")

  notS = (1:n)[-S]
  crit = nout::critWMW(m=m, n=n,alpha=alpha)$crit.vals

  # Ranks of S_Y[i] in (S_X,S_Y[i])
  U_S = sum(S_Y[S]>S_X)
  U_i = sort(sapply(notS, function(i) sum(S_Y[i]>S_X)), decreasing = TRUE)

  # For each k in {1,...,m} consider the worst case scenario
  # (consider the k smallest U_i)
  U.partial = sapply(1:(n-1), function(k) sum(U_i[k:(n-1)]))
  U = c(U.partial, 0)+U_S

  d = ifelse(sum(cumsum(U >= crit) == 1:n)==n, 1, 0)

  return(d)
}



