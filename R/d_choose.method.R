



#' d_choose.method
#'
#' @param S_Y : score vector for the test set
#' @param S_X : score vector for the calibration set
#' @param methods : a character vector indicating one or more of the following local tests to be used in the closed testing procedure:
#' \itemize{
#' \item "simes": Simes' test rejects the global null when
#' \item "storeysimes": Simes' test with Storey estimatore for the proportion of true null hypotheses rejects the global null when ,
#' \item "adjustedfisher": adjusted Fisher test rejects the global null when ,
#' \item "t2": LMPI \eqn{T_2} or Wilcoxon-Mann-Whitney test rejects the global null when ,
#' \item "t3": LMPI \eqn{T_3} test rejects the global null when}
#' Default value is the vector including all possible local tests.
#' @param lambda : parameter involved in the computation of Storey estimator, when method "storeysimes" is used.
#' Default value is set equal to 0.5
#' @param n.exact.f : if \eqn{min\{m,n\}\leq n.exact.f} the critical values of the adjusted Fisher test
#' statistic are computed via permutation. Default value is set equal to 10
#' @param n.exact.t2 : if \eqn{min\{m,n\}\leq n.exact.t2} the critical values of the LMPI \eqn{T_2}
#' statistic are computed via permutation. Default value is set equal to 10
#' @param n.exact.t3 : if \eqn{min\{m,n\}\leq n.exact} the critical values of the LMPI \eqn{T_3}
#' statistic are computed via permutation. Default value is set equal to 10
#' @param B.f : number of iteration to approximate Fisher test distribution via permutation. Default value is set equal to 10^3
#' @param B.t3 : number of iteration to approximate \eqn{T_3} test distribution via permutation. Default value is set equal to 10^3
#' @param seed : seed to ensure reproducible results when approximating the test distribution via permutation
#' @param alpha : significance level. Default value is set equal to 0.1
#'
#'
#' @return A vector of integer numbers corresponding to the \eqn{(1 − \alpha)}-confidence lower bound for the
#' number of true discoveries in closed testing procedure using the local tests chosen by the user.
#' The selection set is trivial, i.e., we are interested in testing all the observations in the test set by default.
#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' d_choose.method(S_Y=Sy, S_X=Sx, methods = "T3")
#' d_choose.method(S_Y=Sy, S_X=Sx, methods = c("simes", "adjustedFisher", "t2", "T3"))
#'
#'
d_choose.method = function(S_X, S_Y, methods = c("adjustedfisher", "t2", "t3", "simes", "storeysimes"), lambda=0.5, n.exact.f=10, n.exact.t2=10, n.exact.t3=10, B.f=10^3, B.t3=10^3, seed=123, alpha=0.1){

  methods = tolower(methods)
  ll = length(methods)
  if(sum(methods %in% c("adjustedfisher", "t2", "t3", "simes", "storeysimes"))!=ll){
    stop("Error: each method should be one of the following options: simes, storeysimes, adjustedfisher, t2 or t3.")
  }

  d.F = NULL; d.S = NULL; d.SS = NULL; d.T2 = NULL; d.T3 = NULL;

  if("adjustedfisher" %in% methods){
    d.F = d_adjustedfisher(S_X=S_X, S_Y=S_Y, n.exact=n.exact.f, B=B.f, seed=seed, alpha=alpha)
  }

  if("simes" %in% methods){
    d.S = d_simes(S_X=S_X, S_Y=S_Y, alpha=alpha)
  }

  if("storeysimes" %in% methods){
    d.SS = d_storeysimes(S_X, S_Y, lambda=lambda, alpha=alpha)
  }

  if("t2" %in% methods){
    d.T2 = d_t2(S_X=S_X, S_Y=S_Y, n.exact=n.exact.t2, alpha=alpha)
  }

  if("t3" %in% methods){
    d.T3 = d_t3(S_X=S_X, S_Y=S_Y,  n.exact=n.exact.t3, B=B.t3, seed=seed, alpha=alpha)
  }

  res = vector()
  if(!is.null(d.S))
    res = c(res,"d.Simes" = d.S)
  if(!is.null(d.SS))
    res = c(res,"d.StoreySimes" = d.SS)
  if(!is.null(d.F))
    res = c(res,"d.Fisher" = d.F)
  if(!is.null(d.T2))
    res = c(res,"d.T2" = d.T2)
  if(!is.null(d.T3))
    res = c(res,"d.T3" = d.T3)

  return(res)
}








