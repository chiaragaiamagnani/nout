


#' crit.T3
#'
#' @description Given the number of observations in the calibration set (\eqn{m}) and the number
#' of observations in the test set (\eqn{n}), it returns the vector of critical
#' values \eqn{U(m,h)} of the Wilcoxon-Mann-Whitney test statistic letting \eqn{h} vary between
#'  \eqn{1} and \eqn{n}, while the size of the first sample is kept fixed equal to \eqn{m}.
#'
#'
#' @param S_Y : score vector for the test set
#' @param S_X : score vector for the calibration set
#' @param alpha : significance level of the local test. Default value is set equal to 0.1
#' @param B : number of interation to approximate \eqn{T_3} distribution via permutation
#' @param seed : seed to ensure reproducible results
#'
#'
#'
#' @return A R object of class *crit.vals.info*, which is a list consisting
#' of the following elements: \itemize{
#' \item m number of test observations
#' \item n number of calibration observations
#' \item alpha significance level of the local test
#' \item crit.val critical value of LMPI \eqn{T_3} statistic.}
#'
#'
#' @export
#'
#'@importFrom foreach %dopar%
#' @examples
#' crits = crit.T3(S_X=runif(10),S_Y=runif(min=0, max=0.5, 5),alpha=0.1)
#'
#'
crit.T3 = function(S_X, S_Y, alpha=0.1, B=10^3, seed = 123){

  set.seed(seed)

  m = length(S_X)
  n = length(S_Y)

  crit.v = foreach::foreach(b = 1:B, .combine=cbind) %dopar% {
    N=m+n
    perm = sample(1:N, m)
    S_Z = c(S_X,S_Y)
    S_Z.perm = c(S_Z[perm], S_Z[-perm])
    U = rank(S_Z.perm)[(m+1):N]-1
    T3 = sum(U^2 + U)
    crit.v = stats::quantile(as.vector(T3), probs = 1-alpha)
    return(crit.v)
  }

  crit = mean(crit.v)

  res = list("m" = m, "n" = n, "crit.val" = crit, "alpha" = alpha)
  class(res) = "crit.vals.info"

  return(res)
}







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
#' @param n.exact : if \eqn{min\{m,n\}\leq n.exact} the critical values of the LMPI \eqn{T_3}
#' statistic are computed via permutation. Default value is set equal to 10.
#' @param B : number of interation to approximate \eqn{T_3} distribution via permutation
#' @param seed : seed to ensure reproducible results
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
d_MannWhitneyk3 = function(S_Y,S_X,n.exact=10,B=10^3, seed=123, alpha=0.1){
  n = length(S_Y)
  m = length(S_X)

  if(min(n,m)>n.exact){
    crit = sapply(1:n, function(h){
      N = m+h
      lambda = m/N
      theta = (4*(2-lambda))/(3*(1-lambda)*lambda)
      variance = 64/(45*N*lambda^3*(1-lambda)^3)
      quantile = stats::qnorm(alpha, mean=(choose(h,2)*choose(m,2))/N*theta+h*(2*h^2-3*h+1)/6+h*(h-1)/2,
                              sd = (choose(h,2)*choose(m,2))/N*sqrt(variance), lower.tail = F)
      return(quantile)
    })
  }
  if(min(n,m)<=n.exact){
    crit = sapply(1:n, function(h) crit.T3(S_X=S_X, S_Y=S_Y[1:h], B=B, alpha=alpha, seed=seed)$crit.val)
  }

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




