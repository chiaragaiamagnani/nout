



#' d_Fisher.corrected
#'
#' @param S_Y : score vector for the test set
#' @param S_X : score vector for the calibration set
#' @param gamma : parameter related to the inflation of the variance of the Fisher test statistic when *p*-values are not independent
#' @param alpha : significance level
#'
#' @return An integer which is the \eqn{(1 − \alpha)}-confidence lower bound for
#' the number of true discoveries in closed testing procedure using
#' the corrected Fisher local test applied to conformal *p*-values, that rejects
#' the global null hypothesis when
#' \deqn{\frac{-2\sum_{i=1}^n\log(p_i)+2n(\sqrt(1+\gamma)-1)}{\sqrt(1+\gamma)}\geq \chi^2_{2n}(1-\alpha),}
#' where \chi^2_{2n}(1-\alpha) is the \eqn{(1-\alpha)} quantile of the chi-squared
#' distribution with \eqn{2n} degree of freedom.
#' The selection set, i.e. the set of hypothesis
#' indices that we are interested in is \eqn{[n]=:\{1,...,n\}} by default.
#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' d_Fisher.corrected(S_Y=Sy, S_X=Sx, alpha=0.1)
d_Fisher.corrected = function(S_Y, S_X, gamma=NULL, alpha = 0.1){

  n = length(S_Y)
  m = length(S_X)

  if(is.null(gamma)){
    gamma = n/m
  }

  if(gamma<0 || gamma==0){
    stop("Error: gamma must be greater than 0.")
  }

  # p-values sorted in decreasing order
  pval = sort(sapply(1:n, function(i) (1+sum(S_X >= S_Y[i]))/(m+1)), decreasing = T)

  Fisher.corrected =
    sapply(1:n,
           function(k){
             (-2*sum(log(pval[1:k]))+2*(n:1)*(sqrt(1+gamma)-1)) / (sqrt(1+gamma))
             })

  # Vector of critical values: the first entry is the critical value referred to
  # the first (up) level in the closed testing (where we test the global null) and
  # the last entry is the critical value referred to the last (bottom) level
  # in the closed testing (where we test the elementary hypotheses).
  # crit = (chi2_{2m}(1-alpha),chi2_{2(m-1)}(1-alpha),..., chi2_{2}(1-alpha))
  crit = stats::qchisq(p=1-alpha, df=2*(n:1))

  # For each k in {1,...,m} consider the worst case scenario, i.e., the scenario
  # in which rejecting the null hypothesis is most difficult.
  # In our case, consider the k largest p-values
  # (the largest p-values are those p-values closer to 1 ->
  # log(sth close to 1) is close to 0 -> the corrected Fisher test statistic is smaller
  # and rejecting the null hypothesis is more difficult)

  d = sum(cumsum(Fisher.corrected >= crit) == 1:n)

  return(d)
}
