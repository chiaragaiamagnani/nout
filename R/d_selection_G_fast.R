
#' d_selection_G2
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
#' X = stats::runif(50)
#' Y = replicate(50, rg2(rnull=runif))
#' res = d_selection_G2(S_Y=Y, S_X=X, B=100)
#' res = d_selection_G2(S_Y=Y, S_X=X, S = c(1:40), g.hat = g2, monotonicity="increasing", B=100)
d_selection_G2 <- function(S_Y, S_X, S=NULL, k=NULL, g.hat=NULL, monotonicity=NULL, prop.F=0.5, alpha=0.1, n_perm=10, B=10^3, B_MC=10^3, seed=123){

  if(!is.null(monotonicity)){
    stopifnot("Error: monotonicity must be either increasing, decreasing"= monotonicity%in%c("decreasing", "increasing"))
  }

  n = as.double(length(S_Y))
  m = as.double(length(S_X))
  N = as.double(m+n)
  s = ifelse(is.null(S), n, length(S))
  S_Z = c(S_X, S_Y)

  if(is.null(g.hat)){
    # Compute the individual statistics for each test point using the input data
    m1 = round(prop.F*m)
    X1=sample(S_X,m1)
    X2=setdiff(S_X,X1)
    g.hat = estimate_g(X1=X1, X2=X2, Y=S_Y, constraint=monotonicity, ker="uniform")
  }

  if(is.null(monotonicity))
    res = d_G_cons2(S_X=S_X, S_Y=S_Y, S=S, g.hat=g.hat, k=k, alpha=alpha, n_perm=n_perm, B=B, seed=seed)
  else{
    decr = ifelse(monotonicity=="increasing", FALSE, TRUE)
    res = d_G_monotone2(S_X=S_X, S_Y=S_Y, S=S, g.hat=g.hat, decr=decr, k=k, alpha=alpha, n_perm=n_perm, B=B, seed=seed)
  }

  return(res)

}



#' d_G_monotone2
#'
#' @param S_Y : test score vector
#' @param S_X :  calibration score vector
#' @param S : selection set in the index test set
#' @param g.hat : it can be either a character ("analytical") or a function denoting the outlier density.
#' If g.hat=="analytical" the test statistics are computed analytically without Monte Carlo estimation.
#' @param decr : logical value indicating whether the outlier distribution is decresing (TRUE)
#' or increasing (FALSE). Default value is FALSE
#' @param k : order of the LMPI test statistic to be specified when g.hat is "analytical"
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
#'
#' @export
#'
#' @examples
#' g2 = function(x, k=2) ifelse(x<1 & x>0, k*x^(k-1), 0)
#' rg2 = function(rnull, k=2) max(rnull(k))
#' m = 10; n=10;
#' X = runif(m)
#' Y = replicate(n, rg2(rnull=runif))
#' res = d_G_monotone2(S_Y=Y, S_X=X, S=c(1:7), decr=FALSE,  g.hat=g2, B=100)
d_G_monotone2 = function(S_X, S_Y, S=NULL, g.hat, decr=F, k=NULL, alpha=0.1, n_perm=10, B=10^3, B_MC = 10^3, seed=123){

  n = length(S_Y)
  m = length(S_X)

  if(is.null(S)){
    # S = [n]
    s = n
    Y.S = sort(S_Y, decreasing = decr)
    ZZ = sapply(length(Y.S):1, function(h) c(S_X, Y.S[1:h]))

  } else {
    s = length(S)
    Y.S = sort(S_Y[S], decreasing = decr)
    notS = setdiff(1:n, S)
    Y.notS = sort(S_Y[notS], decreasing = decr)
    Z.up = sapply(length(Y.notS):1, function(h) c(S_X, Y.S, Y.notS[1:h]))
    Z.down = sapply(length(Y.S):1, function(h) c(S_X, Y.S[1:h]))
    ZZ = c(Z.up, Z.down)
  }

  ## Closed-testing shortcut: sort the test points based on their individual statistics
  ## For each k in {1,...,n} consider the worst-case subset of test points with cardinality k
  tentative.d = 0; l = n
  cont = TRUE
  pval.global = 1
  check = F

  while(cont == TRUE & l>0){
    if(is.character(g.hat)){
      if(g.hat=="analytical"){
        stats_G = sapply(1:(l+m), function(h) (k+1)*k_mom_beta(a=h, b=m+l-h+1, k=k))
      } else{
        cat("Error: g.hat must be either a density function or the string analytical.")
      }

    } else {
      # stats_G = stats_G_j_MC(N=m+l, g=g.hat, B=B_MC)
      stats_G = apply(replicate(B_MC, g.hat(sort(stats::runif(m+l)))) , 1, mean)
    }

    R = stat.G(Z=ZZ[[n-l+1]], m=m, stats_G_vector=stats_G)
    T_wc = sum(R)
    crit = asymptotic.critical.G(m=m, n=l, stats_G_vector=stats_G, alpha=alpha)

    if(l==s){
      T.global = T_wc
      pval.global = compute.global.pvalue(T.obs=T.global, m=m, n=s, local.test="g", stats_G_vector=stats_G,
                                          n_perm=n_perm, B=B, seed=seed)
      check = T
    } else {
      pval.global = pval.global
    }

    if(T_wc >= crit){
      tentative.d = tentative.d+1
      l=l-1
    }
    else{
      tentative.d = tentative.d
      cont=FALSE
    }
  }


  d = ifelse(tentative.d-n+s>0, tentative.d-n+s, 0)

  out = list("lower.bound" = d,
             "global.pvalue" = pval.global,
             "S" = S,
             "selection.p.value" = 1)

}




d_G_monotone_fast = function(S_X, S_Y, S=NULL, g.hat, k=NULL, alpha=0.1, n_perm=10, B=10^3, B_MC = 10^3, seed=123){

  n = length(S_Y)
  m = length(S_X)

  if(is.null(S)){
    # S = [n]
    s = n
    Y.S = sort(S_Y, decreasing = F)
    ZZ = sapply(length(Y.S):1, function(h) c(S_X, Y.S[1:h]))
    R = rank(ZZ)
  } else {
    s = length(S)
    Y.S = sort(S_Y[S], decreasing = F)
    notS = setdiff(1:n, S)
    Y.notS = sort(S_Y[notS], decreasing = F)
    Z.up = sapply(length(Y.notS):1, function(h) c(S_X, Y.S, Y.notS[1:h]))
    Z.down = sapply(length(Y.S):1, function(h) c(S_X, Y.S[1:h]))
    ZZ = c(Z.up, Z.down)
    R = rank(ZZ)
  }

  ## Closed-testing shortcut: sort the test points based on their individual statistics
  ## For each k in {1,...,n} consider the worst-case subset of test points with cardinality k
  tentative.d = 0; pval.global = 1

  for (l in n:1){
    if(is.character(g.hat)){
      if(g.hat=="analytical"){
        stats_G_v = sapply(1:(l+m), function(h) (k+1)*k_mom_beta(a=h, b=m+l-h+1, k=k))
      } else{
        cat("Error: g.hat must be either a density function or the string analytical.")
      }
    } else {
      # stats_G_v = stats_G_j_MC(N=m+l, g=g.hat, B=B_MC)
      stats_G_v = apply(replicate(B, g.hat(sort(stats::runif(m+l)))) , 1, mean)

    }

    T_wc = sum(  stats_G_v[ R[(m+1):(m+l)] ]  )
    crit_l = asymptotic.critical.G(m=m, n=l, stats_G_vector=stats_G_v, alpha=alpha)
    tentative.d = tentative.d  + ( T_wc > crit_l )

    if(l==s){
      T.global = T_wc
      pval.global = compute.global.pvalue(T.obs=T.global, m=m, n=s, local.test="g", stats_G_vector=stats_G_v,
                                          n_perm=n_perm, B=B, seed=seed)
    } else {
      pval.global = pval.global
    }

    if (tentative.d < n - l + 1){ break }
  }

  d = ifelse(tentative.d-n+s>0, tentative.d-n+s, 0)

  out = list("lower.bound" = d,
             "global.pvalue" = pval.global,
             "S" = S,
             "selection.p.value" = 1)

}




#' d_G_cons2
#'
#' @param S_Y : test score vector
#' @param S_X :  calibration score vector
#' @param S : selection set in the index test set
#' @param g.hat : it can be either a character ("analytical") or a function denoting the outlier density.
#' If g.hat=="analytical" the test statistics are computed analytically withuout Monte Carlo estimation.
#' @param k : order of the LMPI test statistic to be specified when g.hat is "analytical"
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
#'
#' @export
#'
#' @examples
#' g2 = function(x, k=2) ifelse(x<1 & x>0, k*x^(k-1), 0)
#' rg2 = function(rnull, k=2) max(rnull(k))
#' m = 10; n=10;
#' X = runif(m)
#' Y = replicate(n, rg2(rnull=runif))
#' res = d_G_cons2(S_Y=Y, S_X=X, g.hat=g2, k=1, B=100)
d_G_cons2 = function(S_X, S_Y, S=NULL, g.hat, k=NULL, alpha=0.1, n_perm=10, B=10^3, B_MC=10^3, seed=123){

  m = as.double(length(S_X))
  n = as.double(length(S_Y))
  N = as.double(n+m)
  s = ifelse(is.null(S), n, length(S))

  tentative.d = 0; l = n
  cont = TRUE
  pval.global = 1

  while(cont == TRUE & l>0){
    if(is.character(g.hat)){
      if(g.hat=="analytical"){
        stats_G = sapply(1:(l+m), function(h) (k+1)*k_mom_beta(a=h, b=m+l-h+1, k=k))
      } else{
        cat("Error: g.hat must be either a density function or the string analytical.")
      }

    } else {
      stats_G = apply(replicate(B_MC, sapply(X=sort(stats::runif(m+l)), FUN=g.hat)), 1, mean)
    }

    if(is.null(S)){

      Rx = sapply(1:n, function(i) rank(c(S_X, S_Y[i]))[m+1]-1)
      range_test_ranks_n = (min(Rx)+1):(max(Rx)+l)
      R = sort(stats_G[range_test_ranks_n], decreasing=F)[1:l]

    } else {
      Rx = sapply(1:n, function(i) rank(c(S_X, S_Y[i]))[m+1]-1)
      Rx.S = Rx[S]
      if(l>=(n-s+1)){
        range_test_ranks_supS = (min(Rx)+1):(max(Rx)+l)
        range_test_ranks_S = range_test_ranks_supS
      } else {
        range_test_ranks_subS = (min(Rx.S)+1):(max(Rx.S)+l)
        range_test_ranks_S = range_test_ranks_subS
      }
      # R = stats_G[range_test_ranks_S]
      R = sort(stats_G[range_test_ranks_n], decreasing = F)[1:l]
    }

    T_wc = sum(R)
    crit = asymptotic.critical.G(m=m, n=l, stats_G_vector=stats_G, alpha=alpha)

    if(l==s){
      T.global = T_wc
      pval.global = compute.global.pvalue(T.obs=T.global, m=m, n=s, local.test="g", stats_G_vector=stats_G,
                                          n_perm=n_perm, B=B, seed=seed)
    } else {
      pval.global = pval.global
    }

    if(T_wc >= crit){
      tentative.d = tentative.d+1
      l=l-1
    } else{
      tentative.d = tentative.d
      cont=FALSE
    }
  }

  d = ifelse(tentative.d-n+s>0, tentative.d-n+s, 0)

  out = list("lower.bound" = d,
             "global.pvalue" = pval.global,
             "S" = S,
             "selection.p.value" = 1)

}



