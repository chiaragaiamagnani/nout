



#' stat.Tk
#'
#' @description It computes the ranks of the test observations in the pooled
#' vector of calibration and test scores. Then, component-wisely compute the power
#' from \eqn{1} to \eqn{k} of the rank vector and component-wisely sum them
#'
#' @param Z : pooled score vector with the first \eqn{m} components corresponding
#' to calibration observations and the last \eqn{n} components corresponding
#' to test observations
#' @param m : calibration sample size
#' @param k : order of the LMPI test statistic
#'
#' @return Given the pooled score vector \eqn{Z=(X,Y)} where \eqn{X} is the calibration score
#' vector and \eqn{Y} is the test score vector, for each observation in the test sample
#' it returns
#' \deqn{T_{k,i}=\sum{j=1}^m\mathbb{1}\{X_j<Y_i\}+\left(\sum{j=1}^m\mathbb{1}\{X_j<Y_i\}\right)^2+\ldots+\left(\sum{j=1}^m\mathbb{1}\{X_j<Y_i\}\right)^k.}
#'
#'
stat.Tk <- function(Z, m, k) {

  N = length(Z)
  n = N-m
  R = rank(Z)[(m+1):N]-1

  if(n==1){
    Tk_i = sum(sapply(1:k, function(h){R^h}))
  }
  else{
    Tk_i = apply(sapply(1:k, function(h){R^h}), MARGIN=1, FUN = sum)
  }

  return(Tk_i)
}



#' asymptotic.critical.Tk
#'
#' @description It computes the \eqn{(1-\alpha)}-quantile of LMPI \eqn{T_k} test statistic
#' based on asymptotic normal approximation.

#' @param m : calibration sample size
#' @param n : test sample size
#' @param k : order of the LMPI test statistic
#' @param correlation_remainder : correlation between the remainder and the $U$-statistic when statistic $T_3$ is used
#' @param alpha : significance level. Default value is set equal to 0.1
#'
#'
#' @return It returns the \eqn{(1-\alpha)}-quantile of LMPI \eqn{T_k} test statistic
#' based on asymptotic normal approximation.
#'
asymptotic.critical.Tk <- function(m, n, k, correlation_remainder, alpha=0.1) {

  stopifnot(k>=1)
  stopifnot((correlation_remainder <=1) & (correlation_remainder >= -1))

  # Compute mean and variance of kernel of Tk
  theta.kernel = compute_theta.Tk(m=m,n=n,k=k)
  variance.kernel = compute_variance.Tk(m=m,n=n,k=k)

  if(k==2 & correlation_remainder==0){
    N = m+n
    # mean remainder
    theta.remainder = ((N-1)*(N+n)/3)#*(choose(m,k)*choose(n,k)/N)
    # variance remainder
    z01_Re = ((m + n)^2*(16*m^4 + 16*n^4 + 30*m^3*(-1 + 2*n) +
                           30*m*n^2*(-1 + 2*n) + m^2*(15 - 60*n + 92*n^2)))/(45*m^4*n^4)
    z10_Re = ((m + n)^2*(16*m^4 + 16*n^4 + 30*m^3*(-1 + 2*n) +
                           30*m*n^2*(-1 + 2*n) + m^2*(75 + 4*n*(-75 + 83*n))))/(45*m^4*n^4)
    variance.remainder = ((choose(n,k) * choose(m,k) / N)^2) * (4/N*(z10_Re*N/m+z01_Re*N/n))

    # theta and variance cosidering the remainder
    theta = theta.kernel+theta.remainder
    variance = variance.kernel+variance.remainder

  } else if (k==2 & correlation_remainder!=0){
    N = m+n
    # mean remainder
    theta.remainder = ((N-1)*(N+n)/3)#*(choose(m,k)*choose(n,k)/N)
    # variance remainder
    z01_Re = ((m + n)^2*(16*m^4 + 16*n^4 + 30*m^3*(-1 + 2*n) +
                           30*m*n^2*(-1 + 2*n) + m^2*(15 - 60*n + 92*n^2)))/(45*m^4*n^4)
    z10_Re = ((m + n)^2*(16*m^4 + 16*n^4 + 30*m^3*(-1 + 2*n) +
                           30*m*n^2*(-1 + 2*n) + m^2*(75 + 4*n*(-75 + 83*n))))/(45*m^4*n^4)
    variance.remainder = ((choose(n,k) * choose(m,k) / N)^2) * (4/N*(z10_Re*N/m+z01_Re*N/n))

    # theta and variance cosidering the remainder
    theta = theta.kernel + theta.remainder
    variance = variance.kernel + variance.remainder +
      2*correlation_remainder*sqrt(variance.kernel*variance.remainder)

  } else if (k>2 || k==1) {
    theta = theta.kernel
    variance = variance.kernel
  }

  critical.value = stats::qnorm(alpha, mean=theta, sd = sqrt(variance), lower.tail = F)
  return(critical.value)
}

#' asymptotic.pvalue.Tk
#'
#' @description It computes the approximated *p*-value of the LMPI \eqn{T_k}
#' test statistic based on the asymptotic normal approximation.
#'
#' @param m : calibration sample size
#' @param n : test sample size
#' @param k : order of the LMPI test statistic
#' @param T.obs : observed value of the test statistic
#' @param correlation_remainder : correlation between the remainder and the $U$-statistic when statistic $T_3$ is used
#'
#'
#' @return It returns the approximated *p*-value of the LMPI \eqn{T_k}
#' test statistic based on the asymptotic normal approximation
#'
asymptotic.pvalue.Tk <- function(m, n, k, T.obs, correlation_remainder) {
  stopifnot(k>=1)
  stopifnot((correlation_remainder <=1) & (correlation_remainder >= -1))

  # Compute mean and variance of kernel of Tk
  theta.kernel = compute_theta.Tk(m=m,n=n,k=k)
  variance.kernel = compute_variance.Tk(m=m,n=n,k=k)

  if(k==2 & correlation_remainder==0){
    N = m+n
    # mean remainder
    theta.remainder = ((N-1)*(N+n)/3)#*(choose(m,k)*choose(n,k)/N)
    # variance remainder
    z01_Re = ((m + n)^2*(16*m^4 + 16*n^4 + 30*m^3*(-1 + 2*n) +
                           30*m*n^2*(-1 + 2*n) + m^2*(15 - 60*n + 92*n^2)))/(45*m^4*n^4)
    z10_Re = ((m + n)^2*(16*m^4 + 16*n^4 + 30*m^3*(-1 + 2*n) +
                           30*m*n^2*(-1 + 2*n) + m^2*(75 + 4*n*(-75 + 83*n))))/(45*m^4*n^4)
    variance.remainder = ((choose(n,2) * choose(m,2) / N)^2) * (4/N*(z10_Re*N/m+z01_Re*N/n))

    # theta and variance cosidering the remainder
    theta = theta.kernel+theta.remainder
    variance = variance.kernel+variance.remainder

  } else if (k==2 & correlation_remainder!=0){

    N = m+n
    # mean remainder
    theta.remainder = ((N-1)*(N+n)/3)#*(choose(m,k)*choose(n,k)/N)
    # variance remainder
    z01_Re = ((m + n)^2*(16*m^4 + 16*n^4 + 30*m^3*(-1 + 2*n) +
                           30*m*n^2*(-1 + 2*n) + m^2*(15 - 60*n + 92*n^2)))/(45*m^4*n^4)
    z10_Re = ((m + n)^2*(16*m^4 + 16*n^4 + 30*m^3*(-1 + 2*n) +
                           30*m*n^2*(-1 + 2*n) + m^2*(75 + 4*n*(-75 + 83*n))))/(45*m^4*n^4)
    variance.remainder = ((choose(n,2) * choose(m,2) / N)^2) * (4/N*(z10_Re*N/m+z01_Re*N/n))

    # theta and variance cosidering the remainder
    theta = theta.kernel + theta.remainder
    variance = variance.kernel + variance.remainder +
      2*correlation_remainder*sqrt(variance.kernel*variance.remainder)

  } else if (k>2 || k==1) {
    theta = theta.kernel
    variance = variance.kernel
  }

  p.value = stats::pnorm(q=T.obs, mean=theta, sd = sqrt(variance), lower.tail = F)

  return(p.value)
}




#' stat.Fisher
#'
#' @description It computes the inverse of conformal *p*-values for each observations in the test set.
#' Then, they are transformed applying the logarithmic function.
#'
#' @param Z : pooled score vector with the first \eqn{m} components corresponding
#' to the calibration observations and the last \eqn{n} components corresponding
#' to the calibration observations
#' @param m : calibration sample size
#'
#'
#' @return Given the pooled score vector \eqn{Z=(X,Y)} where \eqn{X} is the calibration score
#' vector and \eqn{Y} is the test score vector, for each observation in the test sample
#' it returns the following quantity: \deqn{-log(p_i),}
#' where \eqn{p_i = \frac{1+\sum{j=1}^m\mathbb{1}\{X_j<Y_i\}}{m+1}} is the conformal
#' *p*-value of the \eqn{i}th observation in the test set.
#'
#'
stat.Fisher <- function(Z, m) {
  N = length(Z)
  n = N-m
  S_X = Z[1:m]
  S_Y = Z[(m+1):N]
  pval = sapply(1:n, function(i) (1+sum(S_X >= S_Y[i]))/(m+1))
  R = -2 * log(pval)
  return(R)
}


#' asymptotic.critical.Fisher
#'
#' @description It computes the \eqn{(1-\alpha)}-quantile of the adjusted Fisher test statistic
#' based on asymptotic Chi-squared approximation.

#' @param m : calibration sample size
#' @param n : test sample size
#' @param alpha : significance level. Default value is set equal to 0.1
#'
#'
#' @return It returns the \eqn{(1-\alpha)}-quantile of adjusted Fisher test statistic
#' based on asymptotic Chi-squared approximation
#' \deqn{\frac{-2\sum_{i=1}^n\log(p_i) + 2n\sqrt{\gamma+1}-1}{\sqrt{\gamma+1}} \overset{d}{\sim} \chi^2_{2n},}
#' where \eqn{\gamma\in(0,+\infty)} is chosen equal to \eqn{n/m}, \eqn{p_i} are the conformal *p*-values and
#' \eqn{\chi^2_{2n}} is the Chi-squared distribution with \eqn{2n} degrees of freedom.
#'
#'
asymptotic.critical.Fisher <- function(m, n, alpha=0.1) {
  gamma = n/m
  critical.value = sqrt(1+gamma) * stats::qchisq(p=1-alpha, df=2*n) - 2 * (sqrt(1+gamma)-1) * n
  return(critical.value)
}


#' asymptotic.pvalue.Fisher
#'
#' @description It computes the approximated *p*-value of the adjusted Fisher
#' test statistic based on the asymptotic Chi-squared approximation.
#'
#' @param m : calibration sample size
#' @param n : test sample size
#' @param T.obs : observed value of the test statistic
#'
#'
#' @return It returns the approximated *p*-value of the adjusted Fisher
#' test statistic based on the asymptotic Chi-squared approximation
#' \deqn{\frac{-2\sum_{i=1}^n\log(p_i) + 2n\sqrt{\gamma+1}-1}{\sqrt{\gamma+1}} \overset{d}{\sim} \chi^2_{2n},}
#' where \eqn{\gamma\in(0,+\infty)} is chosen equal to \eqn{n/m}, \eqn{p_i} are the conformal *p*-values and
#' \eqn{\chi^2_{2n}} is the Chi-squared distribution with \eqn{2n} degrees of freedom.
#'
#'
asymptotic.pvalue.Fisher <- function(m, n, T.obs) {
  gamma = n/m
  T.obs_shifted = (T.obs+2 * (sqrt(1+gamma)-1) * n)/ sqrt(1+gamma)
  p.value = stats::pchisq(q=T.obs_shifted, df=2*n, lower.tail = F)
  return(p.value)
}




#' perm.crit.T
#'
#' @description It computes permutation \eqn{(1-\alpha)}-quantile of a chosen test statistic.
#'
#' @param m : calibration sample size
#' @param n : test sample size
#' @param k : order of the LMPI test statistic. If \code{NULL} it refers to Fisher test statistic
#' @param alpha : significance level. Default value is set equal to 0.1
#' @param B : number of permutations
#' @param seed : seed to ensure reproducible results
#'
#'
#' @return It returns the \eqn{(1-\alpha)}-quantile of the chosen test statistic obtained via permutation.
#'
#'
perm.crit.T <- function(m, n, k=NULL, alpha=0.1, B=10^3, seed=123){
  set.seed(seed)

  T.v = foreach::foreach(b = 1:B, .combine=cbind) %dopar% {
    N=m+n
    Z = stats::runif(N)
    if(is.null(k)){
      T = sum(stat.Fisher(Z=Z, m=m))
    }
    else {
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
#' @param m : calibration sample size
#' @param n : test sample size
#' @param k : order of the LMPI test statistic. If \code{NULL} it refers to Fisher test statistic
#' @param B : number of permutations
#' @param seed : seed to ensure reproducible results
#'
#'
#' @return It returns the *p*-value for the global null obtained via permutation using the chosen test statistic.
#'
#'
compute.perm.pval <- function(T.obs, m, n, k=NULL, B=10^3, seed=123) {

  set.seed(seed)

  T.v = foreach::foreach(b = 1:B, .combine=cbind) %dopar% {
    N=m+n
    Z = stats::runif(N)
    if(is.null(k))
      T = sum(stat.Fisher(Z=Z, m=m))
    else
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
#' @param k : order of the LMPI test statistic. If \code{NULL} it refers to Fisher test statistic
#' @param n_perm : if \eqn{min(m,n)\leq n_perm} the *p*-value for the global null will be computed via permutation. Default value is 10
#' @param correlation_remainder : correlation between the remainder and the $U$-statistic when statistic $T_3$ is used
#' @param B : number of permutations
#' @param seed : seed to ensure reproducible results
#'
#'
#' @return It returns the *p*-value for the global null, according to the chosen test statistic.
#' The *p*-value is computed via permutation if either the calibration sample size or the test sample size
#' is smaller than \code{n_perm}. Otherwise, it is computed using the asymptotic distribution.
#'
#'
compute.global.pvalue <- function(T.obs, m, n, k=NULL, correlation_remainder, n_perm=10, B=100, seed=321) {

  stopifnot((correlation_remainder <=1) & (correlation_remainder >= -1))

  # permutation p-value for the global null if the sample size is small
  if(min(m,n)<=n_perm){
    if(is.null(k))
      pval.perm = compute.perm.pval(T.obs=T.obs, m=m, n=n, k=k, B=B, seed=seed)
    else
      pval.perm = compute.perm.pval(T.obs=T.obs, m=m, n=n, k=k, B=B, seed=seed)
  }
  # Otherwise, use the asymptotic approximation to compute an approximate p-value for the global null
  else {
    if(is.null(k))
      pval.perm = asymptotic.pvalue.Fisher(m=m, n=n, T.obs=T.obs)
    else
      pval.perm = asymptotic.pvalue.Tk(m=m, n=n, k=k, T.obs=T.obs, correlation_remainder)
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
#' @param alpha : significance level
#' @param k : order of the LMPI test statistic. If \code{NULL} it refers to Fisher test statistic
#' @param correlation_remainder : correlation between the remainder and the $U$-statistic when statistic $T_3$ is used
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
compute.critical.values <- function(m, n, alpha, k=NULL, correlation_remainder, n_perm=10, B=10^3, critical_values=NULL, seed=123){

  stopifnot((correlation_remainder <=1) & (correlation_remainder >= -1))

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
        cat(sprintf("Running permutations...\n"))

        critical.value = perm.crit.T(m=m, n=h, k=k, alpha=alpha, B=B, seed=seed)
      }
    }
    # For large values of m or n compute critical values using the asymptotic approximation of the test statistic
    else {
      if(is.null(k)){
        critical.value = asymptotic.critical.Fisher(m=m, n=h, alpha=alpha)
      } else {
        critical.value = asymptotic.critical.Tk(m=m, n=h, k=k, alpha=alpha, correlation_remainder=correlation_remainder)
      }
    }
    return(critical.value)
  })

  return(crit)
}








#' d_t
#'
#' @description It returns the lower bound for the number of true discoveries in the whole test set
#' obtained with closed testing procedure using LMPI \eqn{T_k} or Fisher local test.
#' No selection in the index test set is performed and the lower bound is computed
#' considering all the observations in the test set.
#'
#' @param S_Y : test score vector
#' @param S_X : calibration score vector
#' @param statistic : parameter indicating the local test to be used in closed testing procedure.
#' It can be either \eqn{T_k} or adjusted Fisher test
#' @param alpha : significance level
#' @param n_perm : if \eqn{min(m,n)\leq n_{perm}} critical values will be computed via permutation. Default value is 10
#' @param B : number of permutation to compute critical values. Default value is 10^3
#' @param critical_values : if not \code{NULL}, a vector of precomputed critical values obtained using
#' the permutation distribution of the test statistic
#' @param correlation_remainder : correlation between the remainder and the $U$-statistic when statistic $T_3$ is used
#' @param seed : seed to ensure reproducible results
#'
#' @return
#' @return A list:
#' \itemize{
#' \item \code{lower_bound}: an integer which is the \eqn{(1 − \alpha)}-confidence lower bound for
#' the number of true discoveries in closed testing procedure using the chosen local test
#' \item \code{global.pvalue}: the global *p*-value, i.e., the *p*-value that closed testing procedure uses to reject the global null
#' }
#'
#' @importFrom foreach %dopar%
#'
#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' d_t(S_Y=Sy, S_X=Sx, statistic="T2", correlation_remainder = 0, alpha=0.1)
#' d_selection_t(S_Y=Sy, S_X=Sx, statistic="T2", alpha=0.1)
#' d_t(S_Y=Sy, S_X=Sx, statistic="T2", correlation_remainder = 1, alpha=0.1)
#'
#' d_t(S_Y=Sy, S_X=Sx, statistic="T3", alpha=0.1)
#' d_selection_t(S_Y=Sy, S_X=Sx, statistic="T3", alpha=0.1)
#'
#' d_t(S_Y=Sy, S_X=Sx, statistic="fisher", alpha=0.1)
#' d_selection_t(S_Y=Sy, S_X=Sx, statistic="fisher", alpha=0.1)
#'
#' d_t(S_Y=Sy, S_X=Sx, statistic="T5", alpha=0.1)
#' d_selection_t(S_Y=Sy, S_X=Sx, statistic="T5", alpha=0.1)
#'
#'
d_t <- function(S_Y, S_X, statistic="T2", alpha=0.1, n_perm=10, B=10^3, critical_values=NULL, correlation_remainder = 0, seed=123){

  statistic = tolower(statistic)

  if(statistic !="fisher"){
    stopifnot(nchar(statistic)>=2)
    cond1 = substring(statistic, 1, 1)=="t"
    cond2 = regmatches(statistic, gregexpr("[[:digit:]]+", statistic))==substring(statistic, 2, nchar(statistic))
    stopifnot(cond1 & cond2)
  }

  if (statistic=="fisher"){
    k=NULL
  } else {
    k = (as.numeric(regmatches(statistic, gregexpr("[[:digit:]]+", statistic))))-1
    stopifnot(k>=1)
  }

  stopifnot((correlation_remainder <=1) & (correlation_remainder >= -1))

  n = length(S_Y)
  m = length(S_X)

  # Compute all critical values for (m,k) from k in {1,...,n}
  crit = compute.critical.values(m=m, n=n, alpha=alpha, k=k, n_perm=n_perm, B=B,
                                 critical_values=critical_values,
                                 correlation_remainder=correlation_remainder, seed=seed)

  # Compute the individual statistics for each test point using the input data
  S_Z = c(S_X, S_Y)
  if(is.null(k)){
    R = stat.Fisher(Z=S_Z, m=m)
  } else{
    R = stat.Tk(Z=S_Z, k=k, m=m)
  }

  ## Closed-testing shortcut: sort the test points based on their individual statistics
  ## For each k in {1,...,n} consider the worst-case subset of test points with cardinality k
  T_i_sorted = sort(R, decreasing = FALSE)
  T_wc = sapply(1:n, function(h) sum(T_i_sorted[1:h]))

  ## Compare the worst-case statistics to the critical values for k in {n,...,1}, starting from the max cardinality
  d = sum(cumsum(rev(T_wc) >= rev(crit)) == 1:n)

  ## Compute p-value for the global null
  T.global = sum(R)
  pval.global = compute.global.pvalue(T.obs=T.global, m=m, n=n, k=k,
                                      correlation_remainder=correlation_remainder,
                                      n_perm=n_perm, B=B, seed=seed)


  ## cat(sprintf("d=%d, p.global=%.3f.\n", d, pval.global))
  ## if( (pval.global > alpha) && (d>0)) {
  ##     cat(sprintf("Inconsistency! d=%d, p.global=%.3f.\n", d, pval.global))
  ## }
  ## if( (pval.global < alpha) && (d==0)) {
  ##     cat(sprintf("Inconsistency! d=%d, p.global=%.3f.\n", d, pval.global))
  ## }

  out = list("lower.bound" = d, "global.pvalue" = pval.global)

  return(out)
}




