



#' d_selection_G
#'
#' @param S_Y : test score vector
#' @param S_X :  calibration score vector
#' @param S : selection set in the index test set
#' @param k : order of the generalized Wilcoxon rank sum test. Classic Wilcoxon test corresponds to \eqn{k=1}
#' @param g.hat : it can be either a character ("analytical") or a function denoting the outlier density.
#' If g.hat=="analytical" the test statistics are computed analytically withuout Monte Carlo estimation.
#' If NULL it is estimated from the data
#' @param monotonicity : character indicating if the outlier density function is monotone increasing or decreasing or neither. Default value is NULL
#' @param prop.F  : proportion of inliers used to estimate the inliser distribution while estimating the outlier density.
#' Default value is 0.5
#' @param alpha : significance level
#' @param n_perm : minimum test sample size needed to use the asymptotic distribution of the test statistic when
#' local.test is either "higher" or "fisher"
#' @param B : number of replications to compute critical values and global *p*-value. Default value is 10^3
#' @param B_MC : number of replications to compute the Shiraishi test statistic
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
#' res = d_selection_G(S_Y=Y, S_X=X, S = c(1:7), g.hat = g2, monotonicity="increasing", B=100)
d_selection_G <- function(S_Y, S_X, S=NULL, k=NULL, g.hat=NULL, monotonicity=NULL, prop.F=0.5, alpha=0.1, n_perm=10, B=10^3, B_MC=10^3, seed=123){

  if(!is.null(monotonicity))
    stopifnot("Error: monotonicity must be either increasing, decreasing"= monotonicity%in%c("decreasing", "increasing"))

  n = as.double(length(S_Y))
  m = as.double(length(S_X))
  N = as.double(m+n)
  s = ifelse(is.null(S), n, length(S))
  S_Z = c(S_X, S_Y)

  if(is.null(g.hat)){

    # Compute the individual statistics for each test point using the input data
    m1 = round(prop.F*m)
    X1 = sample(S_X,m1)
    X2 = setdiff(S_X,X1)
    g.hat = estimate_g(X1=X1, X2=X2, Y=S_Y, constraint=monotonicity, ker="uniform")
    stats_G = sapply(n:1, function(h) apply(replicate(B, g.hat(sort(stats::runif(m+h)))) , 1, mean))


  } else if(is.character(g.hat)){

    if(g.hat=="analytical"){
      stats_G = sapply(n:1, function(l){
        sapply(1:(l+m), function(h) (k+1)*k_mom_beta(a=h, b=m+l-h+1, k=k))})
    } else
      cat("Error: g.hat must be either a density function or the string analytical.")

  } else {
    stats_G = sapply(n:1, function(h) apply(replicate(B, g.hat(sort(stats::runif(m+h)))) , 1, mean))
  }


  if(is.null(monotonicity))
    res = d_G_cons(S_X=S_X, S_Y=S_Y, S=S, stats_G_vector=stats_G, alpha=alpha, n_perm=n_perm, B=B, seed=seed)
  else
    res = d_G_monotone(S_X=S_X, S_Y=S_Y, S=S, stats_G_vector=stats_G, alpha=alpha, n_perm=n_perm, B=B, seed=seed)

  return(res)

}




#' d_G_cons
#'
#' @param S_Y : test score vector
#' @param S_X :  calibration score vector
#' @param S : selection set in the index test set
#' @param stats_G_vector : list of Shiraishi (1985) test statistics for each closed testing level
#' @param alpha : significance level
#' @param n_perm : minimum test sample size needed to use the asymptotic distribution of the test statistic when
#' local.test is either "higher" or "fisher"
#' @param B : number of replications to compute critical values and global *p*-value. Default value is 10^3
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
#' @export
#'
#' @examples
#' g2 = function(x, k=2) ifelse(x<1 & x>0, k*x^(k-1), 0)
#' rg2 = function(rnull, k=2) max(rnull(k))
#' m = 10; n=10;
#' X = runif(m)
#' Y = replicate(n, rg2(rnull=runif))
#' stats_G = sapply(n:1, function(h) apply(replicate(100, g2(sort(runif(m+h)))) , 1, mean))
#' res = d_G_cons(S_Y=Y, S_X=X, stats_G_vector=stats_G, B=100)
d_G_cons = function(S_X, S_Y, S=NULL, stats_G_vector, alpha=0.1, n_perm=10, B=10^3, seed=123){

  m = as.double(length(S_X))
  n = as.double(length(S_Y))
  N = as.double(n+m)
  s = ifelse(is.null(S), n, length(S))

  if(is.null(S)){

    Rx = sapply(1:n, function(i) rank(c(S_X, S_Y[i]))[m+1]-1)
    range_test_ranks_n = sapply(n:1, function(h) (min(Rx)+1):(max(Rx)+h))
    R = sapply(n:1, function(h) sort(stats_G_vector[[h]][range_test_ranks_n[[h]]], decreasing = F))

  } else {
    Rx = sapply(1:n, function(i) rank(c(S_X, S_Y[i]))[m+1]-1)
    Rx.S = Rx[S]
    range_test_ranks_subS = sapply(s:n, function(h) (min(Rx.S)+1):(max(Rx.S)+n-h+1))
    range_test_ranks_supS = sapply(1:s, function(h) (min(Rx)+1):(max(Rx)+n-h+1))
    range_test_ranks_S = c(range_test_ranks_supS, range_test_ranks_subS)
    R = sapply(n:1, function(h) sort(stats_G_vector[[h]][range_test_ranks_S[[h]]], decreasing = F))
  }

  crit = as.double(sapply(1:n, function(h) asymptotic.critical.G(m=m, n=n-h+1, stats_G_vector[[h]], alpha=alpha)))
  T_wc = sapply(1:length(R), function(h) sum(R[[h]][1:h]))

  ## Compare the worst-case statistics to the critical values for k in {n,...,1}, starting from the max cardinality
  tentative.d = as.double(sum(cumsum(rev(T_wc) >= rev(crit)) == 1:n))-n+s
  d = ifelse(tentative.d>0, tentative.d, 0)

  ## Compute p-value for the global null
  T.global = T_wc[s]

  pval.global = compute.global.pvalue(T.obs=T.global, m=m, n=n, local.test="g", stats_G_vector=stats_G_vector[[n-s+1]],
                                      n_perm=n_perm, B=B, seed=seed)

  ## Compute p-value for the selected null
  ## NOTE: this calculation is missing
  pval.selection = 1

  out = list("lower.bound" = d,
             "global.pvalue" = pval.global,
             "S" = S,
             "selection.p.value" = 1)
}









#' d_G_monotone
#'
#' @param S_Y : test score vector
#' @param S_X :  calibration score vector
#' @param S : selection set in the index test set
#' @param stats_G_vector : list of Shiraishi (1985) test statistics for each closed testing level
#' @param alpha : significance level
#' @param n_perm : minimum test sample size needed to use the asymptotic distribution of the test statistic when
#' local.test is either "higher" or "fisher"
#' @param B : number of replications to compute critical values and global *p*-value. Default value is 10^3
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
#' @export
#'
#' @examples
#' g2 = function(x, k=2) ifelse(x<1 & x>0, k*x^(k-1), 0)
#' rg2 = function(rnull, k=2) max(rnull(k))
#' m = 10; n=10;
#' X = runif(m)
#' Y = replicate(n, rg2(rnull=runif))
#' stats_G =  sapply(n:1, function(h) apply(replicate(100, g2(sort(runif(m+h)))) , 1, mean))
#' res = d_G_monotone(S_Y=Y, S_X=X, S=c(1:8), stats_G_vector=stats_G, B=100)
d_G_monotone = function(S_X, S_Y, S=NULL, stats_G_vector, alpha=0.1, n_perm=10, B=10^3, seed=123){

  n = length(S_Y)
  m = length(S_X)

  if(is.null(S)){
    # S = [n]
    s = n
    Y.S = sort(S_Y, decreasing = F)
    ZZ = sapply(length(Y.S):1, function(h) c(S_X, Y.S[1:h]))

  } else {
    s = length(S)
    Y.S = sort(S_Y[S], decreasing = F)
    notS = setdiff(1:n, S)
    Y.notS = sort(S_Y[notS], decreasing = F)
    Z.up = sapply(length(Y.notS):1, function(h) c(S_X, Y.S, Y.notS[1:h]))
    Z.down = sapply(length(Y.S):1, function(h) c(S_X, Y.S[1:h]))
    ZZ = c(Z.up, Z.down)
  }

  #stats_G_test = sapply(length(ZZ):1, function(h) stat.G(Z=ZZ[[h]], m=m, stats_G_vector=stats_G_vector[[h]]))
  ## Closed-testing shortcut: sort the test points based on their individual statistics
  ## For each k in {1,...,n} consider the worst-case subset of test points with cardinality k
  R = sapply(length(ZZ):1, function(h) stat.G(Z=ZZ[[h]], m=m, stats_G_vector=stats_G_vector[[h]]))


  # Compute all critical values for (m,k) from k in {1,...,n}
  crit = as.double(sapply(1:n, function(h) asymptotic.critical.G (m=m, n=n-h+1, stats_G_vector=stats_G_vector[[h]], alpha=alpha)))

  T_wc = sapply(1:length(R), function(h) sum(R[[h]]))


  ## Compare the worst-case statistics to the critical values for k in {n,...,1}, starting from the max cardinality
  tentative.d = as.double(sum(cumsum(rev(T_wc) >= crit) == 1:n))-n+s
  d = ifelse(tentative.d>0, tentative.d, 0)

  ## Compute p-value for the global null
  T.global = T_wc[s]
  pval.global = compute.global.pvalue(T.obs=T.global, m=m, n=s, local.test="g", stats_G_vector=stats_G_vector[[n-s+1]],
                                      n_perm=n_perm, B=B, seed=seed)
  ## Compute p-value for the selected null
  ## NOTE: this calculation is missing
  pval.selection = 1

  out = list("lower.bound" = d,
             "global.pvalue" = pval.global,
             "S" = S,
             "selection.p.value" = 1)

}










