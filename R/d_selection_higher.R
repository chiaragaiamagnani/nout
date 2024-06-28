

# For Lehmann's alternatives we know that the outlier distribution g is monotone

#' d_selection_higher
#'
#' @param S_X :  calibration score vector
#' @param S_Y : test score vector
#' @param S : selection set in the index test set
#' @param local.test : it can be either "wmw" for Wilcoxon rank sum test or "higher" for higher order Wilcoxon rank sum test
#' @param k : order of the generalized Wilcoxon rank sum test. Classic Wilcoxon test corresponds to \eqn{k=1}
#' @param alpha : significance level
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
#' X = runif(50)
#' Y = replicate(50, rg2(rnull=runif))
#' res = d_selection_higher(S_X=X, S_Y=Y, local.test="WMW", B=100)
#' res = d_selection_higher(S_X=X, S_Y=Y, local.test="higher", k=2, S = c(1:40), B=100)
d_selection_higher = function(S_X, S_Y, S=NULL, local.test="wmw", k=NULL, alpha=0.1, n_perm=10, B=10^3, critical_values=NULL, seed=123){

  local.test=tolower(local.test)
  stopifnot(local.test %in% c("wmw", "higher"))

  if(local.test=="wmw") {
    k=1
  } else { stopifnot(k>1 & k%%1==0) }

  if(k==1) local.test="wmw"

  m = as.double(length(S_X))
  n = as.double(length(S_Y))
  N = as.double(n+m)

  if(is.null(S)){
    # S = [n]
    s = n
    Y.S = sort(S_Y, decreasing = F)
    Z = sapply(length(Y.S):1, function(h) c(S_X, Y.S[1:h]))

  } else {
    s = length(S)
    Y.S = sort(S_Y[S], decreasing = F)
    notS = setdiff(1:n, S)
    Y.notS = sort(S_Y[notS], decreasing = F)
    Z.up = sapply(length(Y.notS):1, function(h) c(S_X, Y.S, Y.notS[1:h]))
    Z.down = sapply(length(Y.S):1, function(h) c(S_X, Y.S[1:h]))
    Z = c(Z.up, Z.down)
  }

  ## Closed-testing shortcut: sort the test points based on their individual statistics
  ## For each k in {1,...,n} consider the worst-case subset of test points with cardinality k
  R = sapply(length(Z):1, function(h) stat.Tk(Z=Z[[h]], m=m, k=k))

  # Compute all critical values for (m,k) from k in {1,...,n}
  crit = as.double(compute.critical.values(m=m, n=n, local.test=local.test,
                                           alpha=alpha, k=k, n_perm=n_perm, B=B, critical_values=critical_values, seed=seed))

  T_wc = sapply(1:length(R), function(h) sum(R[[h]]))


  ## Compare the worst-case statistics to the critical values for k in {n,...,1}, starting from the max cardinality
  tentative.d = as.double(sum(cumsum(rev(T_wc) >= rev(crit)) == 1:n))-n+s
  d = ifelse(tentative.d>0, tentative.d, 0)

  ## Compute p-value for the global null
  T.global = T_wc[s]

  pval.global = compute.global.pvalue(T.obs=T.global, m=m, n=s, local.test="higher", k=k, n_perm=n_perm, B=B, seed=seed)

  ## Compute p-value for the selected null
  ## NOTE: this calculation is missing
  pval.selection = 1

  out = list("lower.bound" = d,
             "global.pvalue" = pval.global,
             "S" = S,
             "selection.p.value" = 1)

}



