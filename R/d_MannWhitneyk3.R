
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
#'
#'
d_MannWhitneyk3 = function(S_Y,S_X,alpha){
  n = length(S_Y)
  m = length(S_X)

  crit = sapply(1:n, function(h){
    N = m+h
    lambda = m/N
    theta = (4*(2-lambda))/(3*(1-lambda)*lambda)
    variance = 64/(45*N*lambda^3*(1-lambda)^3)
    quantile = stats::qnorm(alpha, mean=(choose(h,2)*choose(m,2))/N*theta+h*(2*h^2-3*h+1)/6+h*(h-1)/2,
                            sd = (choose(h,2)*choose(m,2))/N*sqrt(variance), lower.tail = F)
    return(quantile)
  })

  S_Z = c(S_X, S_Y)
  R = rank(S_Z)[(m+1):(m+n)]-1
  R2 = R^2
  T3_i = sort(R+R2, decreasing = T)

  # For each k in {1,...,m} consider the worst case scenario
  # (consider the k smallest U_i)
  T3wc = sapply(1:n, function(h) sum(T3_i[h:n]))

  d = sum(cumsum(T3wc >= sort(crit, decreasing = T)) == 1:n)

  return(d)

}




