

#' findd
#'
#' @param method a string indicating the method to be used.
#' Possible values are "Simes", "StoreySimes", "MannWhitney"
#' @param S_X values of calibration scores
#' @param S_Y values of test scores
#' @param S set of selected indices corresponding to the observations to be tested
#' @param lambda parameter involved in the computation of Storey estimator. Default value is set equal to 0.5
#' @param alpha significance level
#'
#' @return An integer which is the \eqn{(1 − \alpha)}-confidence lower bound for the
#' number of true discoveries in closed testing procedure using the selected method as local
#' test.
#'
#' @export
#'
#' @examples
#'
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' findd(method = "MannWhitney", S=1:length(Sy), S_Y=Sy, S_X=Sx, alpha=0.1)
#' findd(method = "Simes", S=1:length(Sy), S_Y=Sy, S_X=Sx, alpha=0.1)
#' findd(method = "StoreySimes", S=1:10, S_Y=Sy, S_X=Sx, alpha=0.1)
#'
findd = function(method, S_X, S_Y, S, lambda = 0.5, alpha){

  method.list = c("Simes", "StoreySimes", "MannWhitney")

  if(!(method %in% method.list)){
    stop("Error: method should be Simes or StoreySimes or MannWhitney")
  }


  n = length(S_Y)
  s = length(S)

  if(s>n){
    stop(sprintf("Error: S contains more items than the entire index set: 1-%s", n))
    #exit("Error: S contains more items than the entire index set: 1-", n)
  }

  if(max(S)>n){
    stop(sprintf("Maximum value of S not contained in the entire index set: 1-%s", n))
  }

  if(min(S)>n){
    stop(sprintf("Minimum value of S not contained in the entire index set: 1-%s", n))
  }

  if(lambda>=1 || lambda<=0){
    stop("lambda must be a value in (0,1)")
  }

  if(alpha>1 || alpha<0){
    stop("alpha must be a value in [0,1]")
  }


  m = length(S_X)

  if(method == "Simes" & s==n){
    d = d_Simes(S_Y=S_Y, S_X=S_X, alpha)
  }

  if(method == "Simes" & s<n){
    d = dselection_Simes(S_Y=S_Y, S_X=S_X, S=S, alpha=alpha)
  }

  if(method == "StoreySimes" & s==n){
    d = d_StoreySimes(S_Y=S_Y, S_X=S_X, lambda=lambda, alpha=alpha)
  }

  if(method == "StoreySimes" & s<n){
    d = dselection_StoreySimes(S_Y=S_Y, S_X=S_X, S=S, lambda=lambda, alpha=alpha)
  }

  if(method == "MannWhitney" & s==n){
    d = d_MannWhitney(S_Y=S_Y, S_X=S_X, alpha=alpha)
  }

  if(method == "MannWhitney" & s<n){
    d = dselection_MannWhitney(S_Y=S_Y, S_X=S_X, S=S, alpha=alpha)
  }

  return(d)

}







