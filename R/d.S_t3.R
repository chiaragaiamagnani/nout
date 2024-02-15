




#' d.S_t3
#'
#'
#'
#'@description For any aribitrary selection subset of the index test set, the function returns a confidence
#' lower bound for the number of true discoveries provided by closed testing procedure using the LMPI \eqn{T_3} local test.
#'
#' @param S_X : score vector for the calibration set
#' @param S_Y : score vector for the test set
#' @param S : vector of selected indices in the test set
#' @param alpha : significance level. Default level is set equal to 0.1
#' @param n.exact : if \eqn{min\{m,n\}\leq n.exact} the critical values of the LMPI \eqn{T_3}
#' statistic are computed via permutation. Default value is set equal to 10
#' @param B : number of iteration to approximate \eqn{T_3} distribution via permutation
#' @param seed : seed to ensure reproducible results
#'
#'
#' @return An integer which is the \eqn{(1 − \alpha)}-confidence lower bound for
#' the number of true discoveries in closed testing procedure using
#' Wilcoxon-Mann-Whitney local test that rejects the global null hypothesis at level \eqn{\alpha}
#' when the test statistic \deqn{\sum_{i=1}^n\sum_{j=1}^m \mathbb{1}\{X_j<Y_i\}} is greater than
#' the critical value corresponding to the significance level \eqn{\alpha}.
#' The selection set can be any arbitrary subset in the index test set
#' according to which units we want to test.
#'
#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' d.S_t3(S_Y=Sy, S_X=Sx, S=c(1:3, 21), alpha=0.1)
#'
#'
#'
d.S_t3 = function(S_X, S_Y, S, alpha=0.1, n.exact=10, B=10^3, seed=123){

  n = length(S_Y)
  m = length(S_X)

  # statistics = ranks of S_Y[i] in (S_X,S_Y)
  S_Z = c(S_X, S_Y)
  R = rank(S_Z)[(m+1):(m+n)]-1
  g = R+R^2

  # Compute critical values
  if(min(n,m)>n.exact){
    cvs = sapply(1:n, function(h){
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
    cvs = sapply(1:n, function(h) nout::perm.crit.t3(S_X=S_X, S_Y=S_Y[1:h], B=B, alpha=alpha, seed=seed)$crit.val)
  }

  res = sumSome::sumStatsPar(g = g, S = S, alpha = alpha, cvs = cvs)

  return(list("d.S_t3"=res$TD, "S"=S))

}


