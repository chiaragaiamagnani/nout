




#' find_d
#'
#' @param X : calibration score vector
#' @param Y : test score vector
#' @param local.test : local test to be used in the closed testing procedure.
#' It can be either "wmw" for Wilcoxon rank sum test, "higher" for higher order Wilcoxon rank sum tests,
#' "fisher" for Fisher's combination test, "g" for the test by Shiraishi (1985), "simes" for Simes' test or
#' "storey" for Simes' test using Storey's estimator for the proportion of true null hypotheses.
#' @param S : selection set in the index test set
#' @param k : positive integer indicating the order of generalized Wilcoxon rank sum test. Default value is NULL
#' @param monotonicity : logical value indicating if the outlier density function is monotone increasing or decreasing or neither.
#' It must be specified when the local test is "g"
#' @param g.hat : it denotes the outlier density. If NULL it is estimated from the data
#' @param alpha : significance level
#' @param prop.F : proportion of inliers used to estimate the inliser distribution while estimating the outlier density.
#' Default value is 0.5
#' @param lambda : parameter to be specified when computing Storey's estimator. Default value is 0.5
#' @param n_perm : minimum test sample size needed to use the asymptotic distribution of the test statistic when
#' local.test is either "higher" or "fisher"
#' @param B : number of replications to compute critical values and global *p*-value. Default value is 10^3
#' @param B_MC : number of replications to compute the Shiraishi test statistic
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
#'
#' @importFrom foreach %dopar%
#'
#' @export
#'
#' @examples
#' g2 = function(x, k=2) ifelse(x<1 & x>0, k*x^(k-1), 0)
#' rg2 = function(rnull, k=2) max(rnull(k))
#'
#' X = runif(10)
#' Y = replicate(10, rg2(rnull=runif))
#' res = find_d(X, Y, local.test="higher", k=3, B=100)
#' res = find_d(X, Y, local.test="g", g.hat = g2, monotonicity="increasing", B=100)
find_d = function(X, Y, local.test = "wmw", S=NULL, k=NULL, monotonicity=NULL, g.hat=NULL, alpha=0.1, prop.F=0.5, lambda=0.5, n_perm=0, B=10^3, B_MC=10^3, critical_values=NULL, seed=123){

  local.test = tolower(local.test)
  stopifnot(local.test %in% c("wmw", "higher", "fisher", "g", "simes", "storey"))

  if(local.test=="higher"){
    stopifnot(k%%1==0 & k>0)
    if(k==1) local.test = "wmw"
  }

  if(local.test=="wmw") k=1


  if(local.test=="wmw" || local.test=="higher"){

    res = d_selection_higher(S_Y=Y, S_X=X, S=S, k=k, alpha=alpha, n_perm=n_perm, B=B, critical_values=critical_values, seed=seed )

  } else if(local.test=="fisher"){

    res = d_selection_fisher(S_Y=Y, S_X=X, S=S, alpha=alpha, n_perm=n_perm, B=B, critical_values=critical_values, seed=seed )

  } else if(local.test=="g"){

    res = d_selection_G2(S_Y=Y, S_X=X, S=S, g.hat=g.hat, monotonicity=monotonicity, prop.F=prop.F, alpha=alpha, n_perm=n_perm, B=B, B_MC=B_MC, seed=seed)

  } else if(local.test=="simes"){

    res = d_selection_simes(S_Y=Y, S_X=X, S=S, alpha=alpha)

  } else if(local.test=="storey"){

    res = d_selection_storey(S_Y=Y, S_X=X, S=S, alpha=alpha, lambda=lambda)

  }

  return(res)

}



