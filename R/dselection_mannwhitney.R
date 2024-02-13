


#' dselection_MannWhitney
#'
#' @param S_Y score vector for the test set
#' @param S_X score vector for the calibration set
#' @param S set of the selected indices
#' @param alpha significance level
#'
#'
#' @return An integer which is the \eqn{(1 − \alpha)}-confidence lower bound for
#' the number of true discoveries in closed testing procedure using
#' Wilcoxon-Mann-Whitney local test applied to conformal *p*-values.
#' The selection set, i.e., the set of hypothesis
#' indices that we are interested in can be chosen.
#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' dselection_MannWhitney(S_Y=Sy, S_X=Sx, alpha = 0.1, S=1:10)
#'
#'
#'
#'
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @useDynLib nout, .registration=TRUE
#'
dselection_MannWhitney = function(S_Y, S_X, S, alpha){

  n = length(S_Y)
  m = length(S_X)
  s = length(S)
  if(s<1)
    stop("Error: selected index set is empty.")

  if(s>n)
    stop("Error: selected index set should be contained in vector 1:length(S_Y).")

  if(s==n)
    d = d_MannWhitney(S_Y = S_Y, S_X = S_X, alpha = alpha)

  if(1<=s & s<n){
    notS = (1:n)[-S]
    crit = sort(nout::crit.WMW(m=m, n=s,alpha=alpha)$crit.vals)

    u_j = sapply(S, function(j) sum(S_Y[j]>S_X))
    v_j = sapply(notS, function(j) sum(S_Y[j]>S_X))

    u0 = max(c(u_j,v_j))
    v0 = min(c(u_j,v_j))
    myu = c(v0, sort(u_j), rep(u0+1, n-s+1))
    myv = c(v0, sort(v_j), rep(u0+1, s+1))

    d = findDiscSum(s=s, m=n, u=myu, v=myv, cs=crit)

  }

  return(d)
}



