#' d_simes
#'
#' @description It returns the lower bound for the number of true discoveries in closed testing procedure
#' using Simes local test.
#'
#' @param S_Y : score vector of test observations
#' @param S_X : score vector of calibration observations
#' @param storey : boolean parameter indicating whether or not using Storey's estimator for the number of true null hypotheses
#' @param alpha : significance level of the local test. Default value is set equal to 0.1.
#'
#' @return An integer which is the \eqn{(1 − \alpha)}-confidence lower bound for the
#' number of true discoveries in closed testing procedure using Simes local
#' test applied to conformal *p*-values.
#' No selection in the index test set is performed and the lower bound is computed
#' considering all the observations in the test set.
#'
#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' d_simes(S_Y=Sy, S_X=Sx)
#' d_simes(S_Y=Sy, storey=T, S_X=Sx)
#'
#'
d_simes = function(S_Y, S_X, storey = F, lambda_storey=0.5, alpha = 0.1){
  n = length(S_Y)
  m = length(S_X)

  pval = sapply(1:n, function(i) (1+sum(S_X >= S_Y[i]))/(m+1))

  if(storey){
    d = hommel::discoveries(hommel::hommel(pval), alpha = alpha)    
  }
  else{
    pval = sort(pval, decreasing = F)
    simes.pval = sapply(1:m, function(i)
      min(pval[i:m]/seq(from=1, to=m-i+1, by=1)))

    # Building the levels of Simes test with Storey estimator
    # pi.not = sapply(1:m, function(i)
    #   (1+sum(pval[i:m]>lambda_storey))/((m-i+1)*(1-lambda_storey)))

    # Building the levels of Simes test with Storey estimator.
    # Storey estimator will be used in the closed testing procedure
    # in every levels except for lowest ones, when the set of
    # considered pvalues has cardinality less than or equal to 2.

    pi.not.highlevels = sapply(1:(m-2), function(i)
      (1+sum(pval[i:m]>lambda_storey))/((m-i+1)*(1-lambda_storey)))
    pi.not = c(pi.not.highlevels,1,1)
    coeff = seq(from = m, to = 1, by = -1)
    thr = alpha/(coeff*pi.not)
    d = sum(cumsum(simes.pval <= thr) == 1:m)
  }

  return(d)
}


