
#' d_MannWhitneyk3
#'
#' @description It returns the lower bound for the number of true discoveries in closed testing procedure
#' using Wilcoxon-Mann-Whitney local test.
#'
#' @param S_Y : score vector for the test set
#' @param S_X : score vector for the calibration set
#' @param alpha : significance level
#'
#' @return An integer which is the \eqn{(1 − \alpha)}-confidence lower bound for
#' the number of true discoveries in closed testing procedure using
#' the LMPI local test corresponding to \eqn{k=3} applied to conformal *p*-values.
#' The selection set, i.e. the set of hypothesis
#' indices that we are interested in is \eqn{[m]=:\{1,...,m\}} by default.
#'
#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' d_MannWhitneyk3(S_Y=Sy, S_X=Sx, alpha=0.1)
#'
#'
d_MannWhitneyk3 = function(S_Y,S_X,alpha){
  n = length(S_Y)
  m = length(S_X)
  k = 3

  crit = sapply(1:n, function(k){
    theta = m/(3*k)*(m-2)
    z11 = 4/45*m^4-4/5*m^3+29/90*m^2-121/90*m-1/5
    z12 = (7*m^4+8*m^3+2*m^2-2*m)/15
    variance = z11+z12/k
    res = stats::qnorm(alpha, mean=theta, sd = sqrt(variance/n),
                                              lower.tail = F)
    return(res)
    })

  U_i = sort(sapply(1:n, function(i) sum(S_Y[i]>S_X)),
             decreasing = TRUE)
  U2_i = U_i^2

  # For each k in {1,...,m} consider the worst case scenario
  # (consider the k smallest U_i)
  Uk3 = sapply(1:n, function(k) sum(U_i[k:n]+U2_i[k:n]))

  d = sum(cumsum(Uk3 >= crit) == 1:n)

  return(d)

}
