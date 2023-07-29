
#' d_MannWhitneyk3
#'
#' @description It returns the lower bound for the number of true discoveries in closed testing procedure
#' using the test statistic \eqn{\sum_{j=1}^n U_j(U_j+1)}, where \eqn{U_j} is the rank of the
#' \eqn{j}th observation in the test set among the vector of calibration units and \eqn{n} is the
#' number of observations in the test set of the test set.
#'
#' @param S_Y : score vector for the test set
#' @param S_X : score vector for the calibration set
#' @param alpha : significance level
#'
#' @return An integer which is the \eqn{(1 − \alpha)}-confidence lower bound for
#' the number of true discoveries in closed testing procedure using
#' the LMPI local test corresponding to \eqn{k=3} applied to conformal *p*-values.
#' The selection set is the set of indices corresponding to the observations in the test set that
#' we want to test. By default, the selection set is \eqn{[n]=:\{1,...,n\}} and corresponds to the
#' entire test set.
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

  crit = sapply(1:n, function(h){
    theta = m/(3*h)*(m+2)
    z11 = 4/45*m^4-4/5*m^3+29/90*m^2-121/90*m-1/5
    z12 = (7*m^4+8*m^3+2*m^2-2*m)/15
    variance = (z11+z12/h)/sqrt(h)
    quantiles = stats::qnorm(alpha, mean=theta, sd = sqrt(variance), lower.tail = F)
    return(quantiles)
    })

  U_i = sapply(1:n, function(i) sum(S_Y[i]>S_X))
  U2_i = U_i^2
  Uk3 = sort(U_i+U2_i, decreasing = T)

  # For each k in {1,...,m} consider the worst case scenario
  # (consider the k smallest U_i)
  Uk3wc = sapply(1:n, function(h) sum(Uk3[h:n]))

  d = sum(cumsum(Uk3wc >= crit) == 1:n)

  return(d)

}
