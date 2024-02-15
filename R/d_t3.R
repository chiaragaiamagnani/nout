



#' perm.crit.T3
#'
#' @description Given \eqn{m} observations in the calibration set and \eqn{n} observations in the test set,
#' it returns the vector of critical values at level \eqn{\alpha} of the LMPI \eqn{T_3} test statistic at level \eqn{\alpha}.
#'
#'
#' @param S_Y : score vector for the test set
#' @param S_X : score vector for the calibration set
#' @param alpha : significance level of the local test. Default value is set equal to 0.1
#' @param B : number of interations to approximate \eqn{T_3} distribution via permutation.
#' Default value is set equal to \eqn{10^3}.
#' @param seed : seed to ensure reproducible results
#'
#'
#'
#' @return A R object of class *crit.val.info*, which is a list consisting
#' of the following elements: \itemize{
#' \item m number of calibration observations
#' \item n number of test observations
#' \item alpha significance level of the local test
#' \item crit.val critical value of LMPI \eqn{T_3} statistic for fixed values of \eqn{m,n} and \eqn{\alpha}.}
#'
#'
#' @export
#'
#' @importFrom foreach %dopar%
#' @examples
#' perm.crit.T3(S_X=runif(10),S_Y=runif(min=0, max=0.5, 5),alpha=0.1)
#'
#'
perm.crit.T3 = function(S_X, S_Y, alpha=0.1, B=10^3, seed = 123){

  set.seed(seed)

  m = length(S_X)
  n = length(S_Y)

  T3.v = foreach::foreach(b = 1:B, .combine=cbind) %dopar% {
    N=m+n
    perm = sample(1:N, m)
    S_Z = c(S_X,S_Y)
    S_Z.perm = c(S_Z[perm], S_Z[-perm])
    U = rank(S_Z.perm)[(m+1):N]-1
    T3 = sum(U^2 + U)
    return(T3)
  }

  crit = stats::quantile(as.vector(T3.v), probs = 1-alpha)

  res = list("m" = m, "n" = n, "crit.val" = crit, "alpha" = alpha)
  class(res) = "crit.val.info"

  return(res)
}






#' d_t3
#'
#' @description It returns the lower bound for the number of true discoveries in closed testing procedure
#' using the LMPI \eqn{T_3} test statistic.
#'
#' @param S_Y : score vector for the test set
#' @param S_X : score vector for the calibration set
#' @param alpha : significance level
#' @param n.exact : if \eqn{min\{m,n\}\leq n.exact} the critical values of the LMPI \eqn{T_3}
#' statistic are computed via permutation. Default value is set equal to 10
#' @param B : number of interation to approximate \eqn{T_3} distribution via permutation
#' @param seed : seed to ensure reproducible results
#'
#' @return An integer which is the \eqn{(1 − \alpha)}-confidence lower bound for
#' the number of true discoveries in closed testing procedure using
#' the LMPI \eqn{T_3} local test that rejects the global null hypothesis at level \eqn{\alpha} for large values of
#' \deqn{\sum_{j=1}^n R_j(R_j+1),} where \eqn{U_j} is the rank of the
#' \eqn{j}th observation in the test set among the pooled score vector \eqn{(S_X,S_Y)} and \eqn{n} is the
#' number of observations in the test set. The selection set is trivial, i.e.,
#' we are interested in testing all the observations in the test set by default.
#'
#'
#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' d_t3(S_Y=Sy, S_X=Sx, alpha=0.1)
#'
#'
#'
#'
d_t3 = function(S_Y, S_X, alpha=0.1, n.exact=10, B=10^3, seed=123){
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
    crit = sapply(1:n, function(h) perm.crit.T3(S_X=S_X, S_Y=S_Y[1:h], B=B, alpha=alpha, seed=seed)$crit.val)
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




