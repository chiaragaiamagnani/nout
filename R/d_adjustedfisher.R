


#' perm.crit.adjustedfisher
#'
#' @description Given \eqn{m} observations in the calibration set and \eqn{n} observations in the test set,
#' it returns the vector of critical values at level \eqn{\alpha} of the adjusted Fisher test statistic
#' \deqn{\frac{-2\sum_{i=1}^n \log(p_i)+2n(\sqrt{1+\gamma}-1)}{\sqrt{1+\gamma}}}
#' where \eqn{\gamma\in(0,+\infty)}.
#'
#'
#' @param S_Y : score vector for the test set
#' @param S_X : score vector for the calibration set
#' @param alpha : significance level of the local test. Default value is set equal to 0.1
#' @param B : number of interations to approximate the adjusted Fisher test statistic distribution via permutation.
#' Default value is set equal to \eqn{10^3}.
#' @param seed : seed to ensure reproducible results
#'
#'
#'
#' @return A R object of class *crit.val.info*, which is a list consisting
#' of the following elements: \itemize{
#' \item m number of test observations
#' \item n number of calibration observations
#' \item alpha significance level of the local test
#' \item crit.val critical value of the adjusted Fisher test statistic.}
#'
#'
#' @export
#'
#' @importFrom foreach %dopar%
#' @examples
#' crits = perm.crit.adjustedfisher(S_X=runif(10),S_Y=runif(min=0, max=0.5, 5),alpha=0.1)
#'
#'
perm.crit.adjustedfisher = function(S_X, S_Y, alpha=0.1, B=10^3, seed = 123){

  set.seed(seed)

  m = length(S_X)
  n = length(S_Y)

  adjFisher.v = foreach::foreach(b = 1:B, .combine=cbind) %dopar% {
    N=m+n
    perm = sample(1:N, m)
    S_Z = c(S_X,S_Y)
    S_X.perm = S_Z[perm]
    S_Y.perm = S_Z[-perm]
    pval = sapply(1:n, function(i) (1+sum(S_X.perm >= S_Y.perm[i]))/(m+1))
    adjFisher.v = (-2*sum(log(pval[1:n]))+2*n*(sqrt(1+n/m)-1)) / (sqrt(1+n/m))
    return(adjFisher.v)
  }

  crit = stats::quantile(as.vector(adjFisher.v), probs = 1-alpha)

  res = list("m" = m, "n" = n, "crit.val" = crit, "alpha" = alpha)
  class(res) = "crit.val.info"

  return(res)
}






#' d_adjustedfisher
#'
#' @description It returns a confidence lower bound for the number of true discoveries provided
#' by closed testing procedure using the adjusted Fisher local test applied to conformal *p*-values.
#'
#' @param S_Y : score vector for the test set
#' @param S_X : score vector for the calibration set
#' @param alpha : significance level. Default level is set equal to 0.1
#' @param n.exact : if \eqn{min\{m,n\}\leq n.exact} the critical values of the adjusted Fisher test
#' statistic are computed via permutation. Default value is set equal to 10
#' @param B : number of iterations to approximate the adjusted Fisher statistic distribution via permutation
#' @param seed : seed to ensure reproducible results
#'
#' @return An integer which is the \eqn{(1 − \alpha)}-confidence lower bound for
#' the number of true discoveries in closed testing procedure using
#' the adjusted Fisher local test applied to conformal *p*-values. The test rejects
#' the global null hypothesis at level \eqn{\alpha} when
#' \deqn{\frac{-2\sum_{i=1}^n\log(p_i)+2n(\sqrt{1+\gamma}-1)}{\sqrt{1+\gamma}\geq \chi^2_{2n}(1-\alpha)},}
#' where *n* is the size of the test set and \eqn{\chi^2_{2n}(1-\alpha)} is the \eqn{(1-\alpha)} quantile of the chi-squared
#' distribution with \eqn{2n} degree of freedom.
#' The selection set is trivial, i.e., we are interested in testing all the observations in the test set by default.
#'
#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' d_adjustedfisher(S_Y=Sy, S_X=Sx, alpha=0.1)
#'
#'
d_adjustedfisher = function(S_Y, S_X, alpha = 0.1, n.exact=10, B=10^3, seed=123){

  n = length(S_Y)
  m = length(S_X)

  # p-values sorted in decreasing order
  pval = sort(sapply(1:n, function(i) (1+sum(S_X >= S_Y[i]))/(m+1)), decreasing = T)

  # Compute adjusted test statistic for each level of the closed testing
  Fisher.adj =
    sapply(n:1,
           function(k){
             (-2*sum(log(pval[1:k]))+2*k*(sqrt(1+k/m)-1)) / (sqrt(1+k/m))
             })

  # Compute critical values
  if(min(n,m)>n.exact){

    # Vector of critical values: the first entry is the critical value referred to
    # the first (up) level in the closed testing (where we test the global null) and
    # the last entry is the critical value referred to the last (bottom) level
    # in the closed testing (where we test the elementary hypotheses).
    # crit = (chi2_{2m}(1-alpha),chi2_{2(m-1)}(1-alpha),..., chi2_{2}(1-alpha))
    crit = stats::qchisq(p=1-alpha, df=2*(n:1))
  }
  if(min(n,m)<=n.exact){
    crit = sapply(1:n, function(h) perm.crit.adjustedfisher(S_X=S_X, S_Y=S_Y[1:h], B=B, alpha=alpha, seed=seed)$crit.val)
  }


  # For each k in {1,...,m} consider the worst case scenario, i.e., the scenario
  # in which rejecting the null hypothesis is most difficult.
  # In our case, consider the k largest p-values
  # (the largest p-values are those p-values closer to 1 ->
  # log(sth close to 1) is close to 0 -> the corrected Fisher test statistic is smaller
  # and rejecting the null hypothesis is more difficult)

  d = sum(cumsum(Fisher.adj >= crit) == 1:n)

  return(d)
}
