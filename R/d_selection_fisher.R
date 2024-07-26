



#' d_selection_fisher
#'
#' @param S_Y : test score vector
#' @param S_X :  calibration score vector
#' @param S : selection set in the index test set
#' @param alpha : significance level
#' @param pvalue_only : logical value. If TRUE, only the global test is performed
#' @param n_perm : minimum test sample size needed to use the asymptotic distribution of the test statistic when
#' local.test is either "higher" or "fisher"
#' @param B : number of replications to compute critical values and global *p*-value. Default value is 10^3
#' @param critical_values : if not \code{NULL}, a vector of precomputed critical values obtained using
#' the permutation distribution of the test statistic
#' @param seed : seed to ensure reproducible results
#'
#' @return A list:
#' \itemize{
#' \item \code{lower_bound}: an integer which is the \eqn{(1 − \alpha)}-confidence lower bound for
#' the number of true discoveries in closed testing procedure using the chosen local test
#' \item \code{S}: the selection set, i.e., the selected subset of the test indices
#' \item \code{global.pvalue}: the global *p*-value, i.e., the *p*-value that closed testing procedure uses to reject the global null
#' \item \code{selection.pvalue}: *p*-value for the selected null
#' }
#' @export
#'
#' @examples
#' g2 = function(x, k=2) ifelse(x<1 & x>0, k*x^(k-1), 0)
#' rg2 = function(rnull, k=2) max(rnull(k))
#'
#' X = runif(10)
#' Y = replicate(10, rg2(rnull=runif))
#' res = d_selection_fisher(S_Y=Y, S_X=X, B=100)
#' res = d_selection_fisher(S_Y=Y, S_X=X, S = c(1:8), B=100)
d_selection_fisher = function(S_X, S_Y, S=NULL, alpha=0.1, pvalue_only=FALSE, n_perm=10, B=10^3, critical_values=NULL, seed=123){

  n = as.double(length(S_Y))
  m = as.double(length(S_X))
  N = as.double(m+n)

  S_Z = c(S_X, S_Y)

  # Compute the individual statistics for each test point using the input data
  R = stat.Fisher(Z=S_Z, m=m)

  if(!pvalue_only){
    # Compute all critical values for (m,k) from k in {1,...,n}
    crit = as.double(compute.critical.values(m=m, n=n, local.test="fisher", alpha=alpha, n_perm=n_perm, B=B,
                                             critical_values=critical_values,
                                             seed=seed))

    if(!is.null(S)) S = as.vector(S)

    # Compute lower bound for S
    res = sumSome::sumStatsPar(g = R, S = S, alpha = alpha, cvs = crit)
    lower.bound = res$TD

    ## Compute p-value for the global null
    T.global = sum(R)

    pval.global = compute.global.pvalue(T.obs=T.global, local.test="fisher", m=m, n=n,
                                        n_perm=n_perm, B=B, seed=seed)

    ## Compute p-value for the selected null
    ## NOTE: this calculation is missing
    pval.selection = 1

  } else {

    T.global = sum(R)

    pval.global = compute.global.pvalue(T.obs=T.global, local.test="fisher", m=m, n=n,
                                        n_perm=n_perm, B=B, seed=seed)
    pval.selection = 1
    lower.bound = 0
  }


  out = list("lower.bound" = lower.bound,
             "global.pvalue" = pval.global,
             "S" = S,
             "selection.p.value" = 1)

}
