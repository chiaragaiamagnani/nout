#' d_selection_simes
#'
#' @description It returns the lower bound for the number of true discoveries in closed testing procedure
#' using Simes local test.
#'
#' @param S_Y : score vector of test observations
#' @param S_X : score vector of calibration observations
#' @param S : set of selected indices in the test set on which compute the lower bound for the number of outliers.
#' If S is missing, dthe global null hypothesis is used
#' @param alpha : significance level of the local test. Default value is set equal to 0.1.
#'
#' @return A list:
#' \itemize{
#' \item \code{lower_bound}: an integer which is the \eqn{(1 − \alpha)}-confidence lower bound for
#' the number of true discoveries in the selected set \eqn{S} in closed testing procedure using Simes' local test
#' \item \code{S}: the selection set, i.e., the selected subset of the test indices
#' \item \code{global.pvalue}: the global *p*-value, i.e., the *p*-value that closed testing procedure uses to reject the global null
#' }
#'
#'
#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' d_selection_simes(S_Y=Sy, S_X=Sx, S=3)
#' d_selection_simes(S_Y=Sy, S_X=Sx, S=c(3, 7:13))
#' d_selection_simes(S_Y=Sy, S_X=Sx)
#' d_simes(S_Y=Sy, S_X=Sx)
#'
d_selection_simes = function(S_Y, S_X, S=NULL, alpha = 0.1){
  n = length(S_Y)
  m = length(S_X)
  pval = sapply(1:n, function(i) (1+sum(S_X >= S_Y[i]))/(m+1))
  hom = hommel::hommel(pval)

  # Compute p-value for the global null
  pval.global = hommel::localtest(hom)

  # Compute p-value for the selected null
  pval.selection = hommel::localtest(hom, ix=S)

  if(is.null(S)){

    # Lower bound
    d = hommel::discoveries(hom, alpha = alpha)

  } else {

    # Lower bound
    d = hommel::discoveries(hom, ix=S, alpha = alpha)

  }

 
  out = list("lower.bound" = d, "global.p.value" = pval.global, "S"=S, "selection.p.value" = pval.selection)

  return(out)
}



#' d_selection_storeysimes
#'
#' @description It returns the lower bound for the number of true discoveries in closed testing procedure
#' using Simes local test with Storey estimator for the proportion of true null hypotheses.
#'
#' @param S_Y : score vector of test observations
#' @param S_X : score vector of calibration observations
#' @param S : set of selected indices in the test set on which compute the lower bound for the number of outliers.
#' If S is missing, dthe global null hypothesis is used
#' @param alpha : significance level of the local test. Default value is set equal to 0.1
#' @param lambda : parameter involved in the computation of Storey estimator. Default value is set equal to 0.5
#'
#'
#' @return An integer which is the \eqn{(1 − \alpha)}-confidence lower bound for the
#' number of true discoveries in closed testing procedure using Simes local
#' test with Storey's estimator for the proportion of true null hypotheses applied to conformal *p*-values.
#' The selection set is trivial, i.e., we are interested in testing all the observations in the test set by default.
#' Then, Storey estimator is computed as
#' \deqn{\hat\pi_0 = \frac{1+\sum_{i=1}^n \mathbb{1}\{p_i>\lambda\}}{n(1-\lambda)}}
#' where \eqn{n} is the test sample size, \eqn{\lambda\in(0,1)} is a tuning parameter
#' and \eqn{p_i} is the *p*-value corresponding to the \eqn{i}th hypothesis in the test set.
#'
#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' d_selection_storeysimes(S_Y=Sy, S_X=Sx)
#' d_selection_storeysimes(S_Y=Sy, S_X=Sx, S=3)
#' d_selection_storeysimes(S_Y=Sy, S_X=Sx, S=c(3, 7:13))
#' d_selection_storeysimes(S_Y=Sy, S_X=Sx)
#' d_storeysimes(S_Y=Sy, S_X=Sx)
#'
#'
#'
# d_selection_storeysimes = function(S_Y, S_X, S=NULL, alpha = 0.1, lambda=0.5){
#   n = length(S_Y)
#   m = length(S_X)
#   pval = sort(sapply(1:n, function(i)
#     (1+sum(S_X >= S_Y[i]))/(m+1)), decreasing=FALSE)
#
#   simes.pval = sapply(1:n, function(i)
#     min(pval[i:n]/seq(from=1, to=n-i+1, by=1)))
#
#   if(!is.null(S)){
#     n = length(S)
#   }
#   # Building the levels of Simes test with Storey estimator
#   # pi.not = sapply(1:n, function(i)
#   #   (1+sum(pval[i:n]>lambda))/((n-i+1)*(1-lambda)))
#
#   # Building the levels of Simes test with Storey estimator.
#   # Storey estimator will be used in the closed testing procedure
#   # in every levels except for lowest ones, when the set of
#   # considered pvalues has cardinality less than or equal to 2.
#
#   pi.not.highlevels = sapply(1:(n-2), function(i)
#     (1+sum(pval[i:n]>lambda))/((n-i+1)*(1-lambda)))
#   pi.not = c(pi.not.highlevels,1,1)
#   coeff = seq(from = n, to = 1, by = -1)
#   thr = alpha/(coeff*pi.not)
#
#   d = sum(cumsum(simes.pval <= thr) == 1:n)
#
#   ## Calculate the p-value for the global null
#   pval.global = min(1, simes.pval[1] * pi.not * n)
#
#   return(list("d"=d, "global.p.value" = pval.global, "S"=S, "pi.not"=pi.not[1]))
# }
#
