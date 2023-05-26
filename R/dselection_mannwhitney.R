


#' dselection_mannwhitney
#'
#' @param S_Y score vector for the test set
#' @param S_X score vector for the calibration set
#' @param alpha significance level
#' @param selection set of the selected indices
#' @param nonselection set of the not selected indices
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
#' dselection_mannwhitney(S_Y=Sy, S_X=Sx, alpha = 0.1, selection = 1:10, nonselection=11:30)
#'
#'
#'
#'
dselection_mannwhitney = function(S_Y,S_X,alpha,selection,nonselection){

  m = length(S_X)
  n = length(S_Y)
  selected.S_Y = S_Y[selection]
  nonselected.S_Y = S_Y[nonselection]
  s = length(selected.S_Y)

  crit = nout::critWMW(m=m, n=s,alpha=alpha)$crit.vals

  ordered.selected.S_Y = sort(selected.S_Y, decreasing = FALSE)

  ordered.selection = sort(selection, decreasing = FALSE)
  ordered.nonselection = sort(nonselection, decreasing = FALSE)

  U_j = sapply(1:s, function(j) sum(selected.S_Y[j]>S_X))
  V_j = sapply(1:(n-s), function(j) sum(nonselected.S_Y[j]>S_X))

  v0 = max(U_j[1], V_j[1])
  tail = rep(min(U_j[s], V_j[n-s]), times=(n-s))
  V = c(V_j, tail)
  U = c(U_j, tail)

  k=1; b=-1; Q=0;

  for(a in 1:n){
    if(U[k+b+1]>=V[a-k-b] || a==1){
      Q=Q+U[k+b+1]
      b=b+1
    }
    else{
      Q=Q+V[a-k-b]
    }
    while(k<=min(s,a) & Q>crit[n-a+1]){
      if(b>0){
        b=b-1
      }
      else{
        Q=Q+U[k+1]-V[a-k]
      }
      k=k+1
    }
  }

  return(s-k+1)
}



