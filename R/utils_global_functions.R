


#' perm.crit.T
#'
#' @description It computes permutation \eqn{(1-\alpha)}-quantile of a chosen test statistic.
#'
#' @param m : calibration sample size
#' @param n : test sample size
#' @param local.test : local test to be used in the closed testing procedure. Default value is Wilcoxon test
#' @param k : order of the LMPI test statistic. If \code{NULL} it refers to Fisher test statistic
#' @param stats_G_vector : vector of elementary test statistics to perform the test in Shiraishi (1985). If NULL it will be computed in d_t using B_MC iterations
#' @param alpha : significance level. Default value is set equal to 0.1
#' @param B : number of permutations
#' @param seed : seed to ensure reproducible results
#'
#'
#' @return It returns the \eqn{(1-\alpha)}-quantile of the chosen test statistic obtained via permutation.
#'
#'
perm.crit.T <- function(m, n, local.test="wmw", k=NULL, stats_G_vector=NULL, alpha=0.1, B=10^3, seed=123){
  set.seed(seed)

  local.test = tolower(local.test)
  stopifnot(local.test %in% c("wmw", "higher", "fisher", "g"))

  if(local.test=="higher"){
    stopifnot(k%%1==0 & k>0)
    if(k==1) local.test = "wmw"
  }
  if(local.test=="wmw") k=1

  T.v = foreach::foreach(b = 1:B, .combine=cbind) %dopar% {
    N=m+n
    Z = stats::runif(N)
    if(local.test=="fisher"){
      T = sum(stat.Fisher(Z=Z, m=m))
    } else if(local.test=="g"){
      T = sum(stat.G(Z=Z, m=m, stats_G_vector=stats_G_vector))
    } else if(local.test=="higher" || local.test=="wmw"){
      T = sum(stat.Tk(Z=Z, k=k, m=m))
    }
  }

  # Empirical quantile with finite-sample adjustment
  idx = ceiling(B*(1-alpha)*(1+1/B))
  if(idx <= B) {
    crit = sort(as.vector(T.v))[idx]
  } else {
    crit = Inf
  }

  return(crit)
}



#' compute.perm.pval
#'
#' @description It computes the *p*-value for the global null obtained via permutation using a chosen test statistic.
#'
#' @param T.obs : observed value of the chosen test statistic
#' @param local.test : local test to be used in the closed testing procedure. Default value is Wilcoxon test
#' @param m : calibration sample size
#' @param n : test sample size
#' @param k : order of the LMPI test statistic. If \code{NULL} it refers to Fisher test statistic
#' @param stats_G_vector : vector of elementary test statistics to perform the test in Shiraishi (1985). If NULL it will be computed in d_t using B_MC iterations
#' @param B : number of permutations
#' @param seed : seed to ensure reproducible results
#'
#'
#' @return It returns the *p*-value for the global null obtained via permutation using the chosen test statistic.
#'
#'
compute.perm.pval <- function(T.obs, local.test="wmw", m, n, k=NULL, stats_G_vector=NULL, B=10^3, seed=123) {

  set.seed(seed)

  local.test = tolower(local.test)
  stopifnot(local.test %in% c("wmw", "higher", "fisher", "g"))

  if(local.test=="higher"){
    stopifnot(k%%1==0 & k>0)
    if(k==1) local.test = "wmw"
  }
  if(local.test=="wmw") k=1

  T.v = foreach::foreach(b = 1:B, .combine=cbind) %dopar% {
    N=m+n
    Z = stats::runif(N)
    if(local.test=="fisher")
      T = sum(stat.Fisher(Z=Z, m=m))
    else if(local.test=="g")
      T = sum(stat.G(Z=Z, m=m, stats_G_vector=stats_G_vector))
    else if(local.test=="higher" || local.test=="wmw")
      T = sum(stat.Tk(Z=Z, k=k, m=m))
    return(T)
  }

  # Compute the permutation p-value
  pval = (1+sum(T.v >= T.obs)) / (1 + length(T.v))

  return(pval)
}



#' compute.global.pvalue
#'
#' @description It computes the *p*-value for the global null according to the chosen test statistic
#'
#' @param T.obs : observed value of the chosen test statistic
#' @param m : calibration sample size
#' @param n : test sample size
#' @param local.test : local test to be used in the closed testing procedure. Default value is Wilcoxon test
#' @param k : order of the LMPI test statistic. If \code{NULL} it refers to Fisher test statistic
#' @param stats_G_vector : vector of elementary test statistics to perform the test in Shiraishi (1985). If NULL it will be computed in d_t using B_MC iterations
#' @param n_perm : if \eqn{min(m,n)\leq n_perm} the *p*-value for the global null will be computed via permutation. Default value is 10
#' @param B : number of permutations
#' @param seed : seed to ensure reproducible results
#'
#'
#' @return It returns the *p*-value for the global null, according to the chosen test statistic.
#' The *p*-value is computed via permutation if either the calibration sample size or the test sample size
#' is smaller than \code{n_perm}. Otherwise, it is computed using the asymptotic distribution.
#'
#'
compute.global.pvalue <- function(T.obs, m, n, local.test="wmw", k=NULL, stats_G_vector=NULL, n_perm=10, B=100, seed=321) {

  local.test = tolower(local.test)
  stopifnot(local.test %in% c("wmw", "higher", "fisher", "g"))

  if(local.test=="higher"){
    stopifnot(k%%1==0 & k>0)
    if(k==1) local.test = "wmw"
  }
  if(local.test=="wmw") k=1

  # permutation p-value for the global null if the sample size is small
  if(min(m,n)<=n_perm){
    if(local.test=="fisher")
      pval.perm = compute.perm.pval(T.obs=T.obs, local.test="fisher", m=m, n=n, k=k, B=B, seed=seed)
    else if(local.test=="g")
      pval.perm = compute.perm.pval(T.obs=T.obs, local.test=="g", m=m, n=n, k=k, stats_G_vector=stats_G_vector, B=B, seed=seed)
    else if(local.test=="higher" || local.test=="wmw")
      pval.perm = compute.perm.pval(T.obs=T.obs, local.test="higher", m=m, n=n, k=k, B=B, seed=seed)
  }
  # Otherwise, use the asymptotic approximation to compute an approximate p-value for the global null
  else {
    if(local.test=="fisher")
      pval.perm = asymptotic.pvalue.Fisher(m=m, n=n, T.obs=T.obs)
    else if(local.test=="g")
      pval.perm = asymptotic.pvalue.G(m=m, n=n, T.obs=T.obs, stats_G_vector = stats_G_vector)
    else if(local.test=="higher" || local.test=="wmw")
      pval.perm = asymptotic.pvalue.Tk(m=m, n=n, k=k, T.obs=T.obs)
  }
  return(pval.perm)
}



#' compute.critical.values
#'
#' @description It computes the vector of critical values for a chosen test statistic
#' at significance level \eqn{\alpha}.
#'
#' @param m : calibration size
#' @param n : test size
#' @param local.test : local test to be used in the closed testing procedure. Default value is Wilcoxon test
#' @param alpha : significance level
#' @param k : order of the LMPI test statistic. If \code{NULL} it refers to Fisher test statistic
#' @param stats_G_vector : vector of elementary test statistics to perform the test in Shiraishi (1985). If NULL it will be computed in d_t using B_MC iterations
#' @param n_perm : if \eqn{min(m,n)\leq n_perm} critical values will be computed via permutation. Default value is 10
#' @param B : number of permutation to compute critical values. Default value is 10^3
#' @param critical_values : if not \code{NULL}, a vector of precomputed critical values obtained using
#' the permutation distribution of the test statistic
#' @param seed : seed to ensure reproducible results
#'
#'
#' @return A vector of critical values for a test statistic chosen among LMPI \eqn{T_k} or Fisher statistics
#' at significance level \eqn{\alpha} with calibration size \eqn{m} fixed for each level of closed testing.
#'
#'
compute.critical.values <- function(m, n, local.test="wmw", alpha, k=NULL, stats_G_vector=NULL, n_perm=10, B=10^3, critical_values=NULL, seed=123){

  local.test = tolower(local.test)
  stopifnot(local.test %in% c("wmw", "higher", "fisher", "g"))

  if(local.test=="higher"){
    stopifnot(k%%1==0 & k>0)
    if(k==1) local.test = "wmw"
  }
  if(local.test=="wmw") k=1

  crit = sapply(1:n, function(h) {
    # For small values of m and n compute critical values via permutation
    if(min(m,h)<=n_perm) {
      found.value = FALSE
      # In order to avoid repeating computation of precomputed critical values that are saved in "tables" folder
      if(!is.null(critical_values)) {
        if(length(critical_values)>=h) {
          critical.value = critical_values[h]
          found.value = TRUE
        }
      }
      if(!found.value) {
        #cat(sprintf("Running permutations...\n"))

        critical.value = as.double(perm.crit.T(m=m, n=h, k=k, local.test=local.test, alpha=alpha, B=B, seed=seed))
      }
    }
    # For large values of m or n compute critical values using the asymptotic approximation of the test statistic
    else {
      if(local.test=="fisher"){
        critical.value = as.double(asymptotic.critical.Fisher(m=m, n=h, alpha=alpha))
      } else if(local.test=="g"){
        critical.value = as.double(asymptotic.critical.G(m=m, n=h, stats_G_vector=stats_G_vector, alpha=alpha)) # QUI
      } else if(local.test=="higher" || local.test=="wmw"){
        critical.value = as.double(asymptotic.critical.Tk(m=m, n=h, k=k, alpha=alpha))
      }
    }
    return(as.double(critical.value))
  })
  return(crit)
}




