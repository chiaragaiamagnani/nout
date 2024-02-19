
#' stat.T2
#'
#' @description It computes the ranks of the test observations in the pooled
#' vector of calibration and test scores.
#'
#' @param Z : pooled score vector with the first \eqn{m} components corresponding
#' to the calibration observations and the last \eqn{n} components corresponding
#' to the calibration observations
#' @param m : calibration sample size
#'
#'
#' @return Given the pooled score vector \eqn{Z=(X,Y)} where \eqn{X} is the calibration score
#' vector and \eqn{Y} is the test score vector, for each observation in the test sample
#' it returns its rank:\deqn{\sum{j=1}^m\mathbb{1}\{X_j<Y_i\}.}
#'
#'
stat.T2 <- function(Z, m) {
  N = length(Z)
  n = N-m
  R = rank(Z)[(m+1):N]-1
  return(R)
}


#' asymptotic.critical.T2
#'
#' @description It computes the \eqn{(1-\alpha)}-quantile of LMPI \eqn{T_2} (Wilcoxon-Mann-Whitney) test statistic
#' based on asymptotic normal approximation.

#' @param m : calibration sample size
#' @param n : test sample size
#' @param alpha : significance level. Default value is set equal to 0.1
#'
#'
#' @return It returns the \eqn{(1-\alpha)}-quantile of LMPI \eqn{T_2} test statistic
#' based on asymptotic normal approximation
#' \deqn{N\left(\frac{n*(m+n+1)}{2}, \frac{n*m*(n+m+1)}{12}\right).}
#'
#'
asymptotic.critical.T2 <- function(m, n, alpha=0.1) {
  # to add continuity correction modify the mean as: mean=((n*(m+n+1)-1)/2)
  critical.value = stats::qnorm(alpha, mean=((n*(m+n+1))/2), sd = sqrt(n*m*(n+m+1)/12), lower.tail = F)
  return(critical.value)
}


#' asymptotic.pvalue.T2
#'
#' @description It computes the approximated *p*-value of the LMPI \eqn{T_2} (Wilcoxon-Mann-Whitney)
#' test statistic based on the asymptotic normal approximation.
#'
#' @param m : calibration sample size
#' @param n : test sample size
#' @param T.obs : observed value of the test statistic
#'
#'
#' @return It returns the approximated *p*-value of the LMPI \eqn{T_2}
#' test statistic based on the asymptotic normal approximation
#' \deqn{N\left(\frac{n*(m+n+1)}{2}, \frac{n*m*(n+m+1)}{12}\right).}
#'
#'
asymptotic.pvalue.T2 <- function(m, n, T.obs) {
  # to add continuity correction modify the mean as: mean=((n*(m+n+1)-1)/2)
  p.value = stats::pnorm(q=T.obs, mean=((n*(m+n+1))/2), sd = sqrt(n*m*(n+m+1)/12), lower.tail = F)
  return(p.value)
}



#' stat.T3
#'
#' @description It computes the ranks of the test observations in the pooled
#' vector of calibration and test scores. Then, component-wise sum the vector
#' of ranks to the vector of the squared ranks.
#'
#' @param Z : pooled score vector with the first \eqn{m} components corresponding
#' to the calibration observations and the last \eqn{n} components corresponding
#' to the calibration observations
#' @param m : calibration sample size
#'
#'
#' @return Given the pooled score vector \eqn{Z=(X,Y)} where \eqn{X} is the calibration score
#' vector and \eqn{Y} is the test score vector, for each observation in the test sample
#' it returns the following quantity: \deqn{R_i + (R_i)^2,}
#' where \eqn{R_i=\sum{j=1}^m\mathbb{1}\{X_j<Y_i\}} is the rank of the \eqn{i}th observation in the test set.
#'
#'
stat.T3 <- function(Z, m) {
  N = length(Z)
  n = N-m
  R = rank(Z)[(m+1):N]-1
  R2 = R^2
  return(R2+R)
}


#' asymptotic.critical.T3
#'
#' @description It computes the \eqn{(1-\alpha)}-quantile of LMPI \eqn{T_3} test statistic
#' based on asymptotic normal approximation.

#' @param m : calibration sample size
#' @param n : test sample size
#' @param alpha : significance level. Default value is set equal to 0.1
#'
#'
#' @return It returns the \eqn{(1-\alpha)}-quantile of LMPI \eqn{T_3} test statistic
#' based on asymptotic normal approximation
#' \deqn{N \left(\frac{\choose(n,2)*\choose(m,2)*\\theta}{N}+\frac{n*(2*n^2-3*n+1)}{6}+\frac{n*(n-1)}{2},
#' \frac{\choose(n,2)*\choose(m,2))}{N}*\frac{64}{(45*N*\lambda^3*(1-\lambda)^3)\right)}.}
#'
#'
asymptotic.critical.T3 <- function(m, n, alpha=0.1) {
  N = m+n
  lambda = m/N
  theta = (4*(2-lambda))/(3*(1-lambda)*lambda)
  variance = 64/(45*N*lambda^3*(1-lambda)^3)
  critical.value = stats::qnorm(alpha, mean=(choose(n,2)*choose(m,2))/N*theta+n*(2*n^2-3*n+1)/6+n*(n-1)/2,
                                sd = (choose(n,2)*choose(m,2))/N*sqrt(variance), lower.tail = F)
  return(critical.value)
}


#' asymptotic.pvalue.T3
#'
#' @description It computes the approximated *p*-value of the LMPI \eqn{T_3}
#' test statistic based on the asymptotic normal approximation.
#'
#' @param m : calibration sample size
#' @param n : test sample size
#' @param T.obs : observed value of the test statistic
#'
#'
#' @return It returns the approximated *p*-value of the LMPI \eqn{T_3}
#' test statistic based on the asymptotic normal approximation
#' \deqn{N \bigg(\frac{\choose(n,2)*\choose(m,2)*\theta}{N}+\frac{n*(2*n^2-3*n+1)}{6}+\frac{n*(n-1)}{2},
#' \frac{\choose(n,2)*\choose(m,2))}{N}*\frac{64}{(45*N*\lambda^3*(1-\lambda)^3)\bigg)}.}
#'
#'
asymptotic.pvalue.T3 <- function(m, n, T.obs) {
  N = m+n
  lambda = m/N
  theta = (4*(2-lambda))/(3*(1-lambda)*lambda)
  variance = 64/(45*N*lambda^3*(1-lambda)^3)
  p.value = stats::pnorm(q=T.obs, mean=(choose(n,2)*choose(m,2))/N*theta+n*(2*n^2-3*n+1)/6+n*(n-1)/2,
                         sd = (choose(n,2)*choose(m,2))/N*sqrt(variance), lower.tail = F)
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
#' @param stat.func : test statistic of which compute critical values.
#' Can be either "stat.T2" or "stat.T2" or "stat.Fisher"
#' @param alpha : significance level. Default value is set equal to 0.1
#' @param B : number of permutations
#' @param seed : seed to ensure reproducible results
#'
#'
#' @return It returns the \eqn{(1-\alpha)}-quantile of the \code{stat.func} test statistic obtained via permutation.
#'
#'
perm.crit.T <- function(m, n, stat.func, alpha=0.1, B=10^3, seed=123){
    set.seed(seed)

    T.v = foreach::foreach(b = 1:B, .combine=cbind) %dopar% {
        N=m+n
        Z = stats::runif(N)
        T = sum(stat.func(Z, m))
        return(T)
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
#' @param stat.func : test statistic of which compute critical values.
#' It can be either "stat.T2" or "stat.T2" or "stat.Fisher"
#' @param B : number of permutations
#' @param seed : seed to ensure reproducible results
#'
#' @keywords internal
#'
#' @return It returns the *p*-value for the global null obtained via permutation using \code{stat.func} test statistic.
#'
#'
compute.perm.pval <- function(T.obs, m, n, stat.func, B=10^3, seed=123) {
    set.seed(seed)

    T.v = foreach::foreach(b = 1:B, .combine=cbind) %dopar% {
        N=m+n
        Z = stats::runif(N)
        T = sum(stat.func(Z, m))
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
#' @param stat.func : test statistic of which compute critical values.
#' It can be either "stat.T2" or "stat.T2" or "stat.Fisher"
#' @param asymptotic.pvalue.func : asymptotic distribution to compute the *p*-value for the global null
#' @param n_perm : if \eqn{min(m,n)\leq n_perm} the *p*-value for the global null will be computed via permutation. Default value is 10
#' @param B : number of permutations
#' @param seed : seed to ensure reproducible results
#'
#' @keywords internal
#'
#' @return It returns the *p*-value for the global null, according to the chosen test statistic.
#' The *p*-value is computed via permutation if either the calibration sample size or the test sample size
#' is smaller than \code{n_perm}. Otherwise, it is computed using the asymptotic distribution.
#'
#'
compute.global.pvalue <- function(T.obs, m, n, stat.func, asymptotic.pvalue.func, n_perm=10, B=100, seed=321) {

    # permutation p-value for the global null if the sample size is small
    if(min(m,n)<=n_perm){
        pval.perm = compute.perm.pval(T.obs=T.obs, m=m, n=n, stat.func=stat.func, B=B, seed=seed)
      }
    # Otherwise, use the asymptotic approximation to compute an approximate p-value for the global null
    else {
        pval.perm = asymptotic.pvalue.func(m=m, n=n, T.obs=T.obs)
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
#' @param stat.func : test statistic of which compute critical values
#' @param asymptotic.critical.func : asymptotic distribution of the chosen test statistic
#' @param n_perm : if \eqn{min(m,n)\leq n_perm} critical values will be computed via permutation. Default value is 10
#' @param B : number of permutation to compute critical values. Default value is 10^3
#' @param critical_values : if not \code{NULL}, a vector of precomputed critical values obtained using
#' the permutation distribution of the test statistic
#' @param seed : seed to ensure reproducible results
#'
#' @keywords internal
#'
#' @return A vector of critical values for a test statistic chosen among \eqn{T_2, T_3} or Fisher
#' at significance level \eqn{\alpha} with calibration size \eqn{m} fixed for each level of closed testing.
#'
#'
compute.critical.values <- function(m, n, alpha, stat.func, asymptotic.critical.func, n_perm=10, B=10^3, critical_values=NULL, seed=123){

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
                critical.value = perm.crit.T(m, h, stat.func=stat.func, alpha=alpha, B=B, seed=seed)
            }
        }
      # For large values of m or n compute critical values using the asymptotic approximation of the test statistic
      else {
            critical.value = asymptotic.critical.func(m, h, alpha)
        }
        return(critical.value)
    })

    return(crit)
}





#' d_t
#'
#' @description It returns the lower bound for the number of true discoveries in the whole test set
#' obtained with closed testing procedure using LMPI \eqn{T_2} (Wilcoxon-Mann-Whitney)
#' or LMPI \eqn{T_3} or Fisher local test. No selection in the index test set is performed and the lower bound is computed
#' considering all the observations in the test set.
#'
#' @param S_Y : test score vector
#' @param S_X : calibration score vector
#' @param statistic : parameter indicating the local test to be used in closed testing procedure.
#' It can be either \eqn{T_2, T_3} or adjusted Fisher test
#' @param alpha : significance level
#' @param n_perm : if \eqn{min(m,n)\leq n_{perm}} critical values will be computed via permutation. Default value is 10
#' @param B : number of permutation to compute critical values. Default value is 10^3
#' @param critical_values : if not \code{NULL}, a vector of precomputed critical values obtained using
#' the permutation distribution of the test statistic
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
#' d_t(S_Y=Sy, S_X=Sx, statistic="T2", alpha=0.1)
#' d_t(S_Y=Sy, S_X=Sx, statistic="T3", alpha=0.1)
#' d_t(S_Y=Sy, S_X=Sx, statistic="fisher", alpha=0.1)
#'
#'
d_t <- function(S_Y, S_X, statistic="T2", alpha=0.1, n_perm=10, B=10^3, critical_values=NULL, seed=123){

    stopifnot(statistic %in% c("T2", "T3", "fisher"))

    if(statistic=="T2") {
        stat.func = stat.T2
        asymptotic.critical.func = asymptotic.critical.T2
        asymptotic.pvalue.func = asymptotic.pvalue.T2
    } else if (statistic=="T3") {
        stat.func = stat.T3
        asymptotic.critical.func = asymptotic.critical.T3
        asymptotic.pvalue.func = asymptotic.pvalue.T3
    } else if (statistic=="fisher") {
        stat.func = stat.Fisher
        asymptotic.critical.func = asymptotic.critical.Fisher
        asymptotic.pvalue.func = asymptotic.pvalue.Fisher
    }

    n = length(S_Y)
    m = length(S_X)

    # Compute all critical values for (m,k) from k in {1,...,n}
    crit = compute.critical.values(m, n, alpha, stat.func, asymptotic.critical.func, n_perm=n_perm, B=B, critical_values=critical_values, seed=seed)

    # Compute the individual statistics for each test point using the input data
    S_Z = c(S_X, S_Y)
    R = stat.func(S_Z, m)

    ## Closed-testing shortcut: sort the test points based on their individual statistics
    ## For each k in {1,...,n} consider the worst-case subset of test points with cardinality k
    T_i_sorted = sort(R, decreasing = FALSE)
    T_wc = sapply(1:n, function(h) sum(T_i_sorted[1:h]))

    ## Compare the worst-case statistics to the critical values for k in {n,...,1}, starting from the max cardinality
    d = sum(cumsum(rev(T_wc) >= rev(crit)) == 1:n)

    ## Compute p-value for the global null
    T.global = sum(R)
    pval.global = compute.global.pvalue(T.obs=T.global, m=m, n=n, stat.func=stat.func,
                                        asymptotic.pvalue.func=asymptotic.pvalue.func, n_perm=n_perm, B=B, seed=seed)


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
