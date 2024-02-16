
# Returns vector of values of the LMPI T2 statistic for each test observation
stat.T2 <- function(Z, m) {
  N = length(Z)
  n = N-m
  R = rank(Z)[(m+1):N]-1
  return(R)
}

# Returns the approximated (1-alpha)-quantile of LMPI T2 based on asymptotic normal approximation
asymptotic.critical.T2 <- function(m, n, alpha) {
  # to add continuity correction modify the mean as: mean=((n*(m+n+1)-1)/2)
  critical.value = stats::qnorm(alpha, mean=((n*(m+n+1))/2), sd = sqrt(n*m*(n+m+1)/12), lower.tail = F)
  return(critical.value)
}

# Returns the approximated pvalue of theLMPI T2 statistic based on asymptotic normal approximation
asymptotic.pvalue.T2 <- function(m, n, T.obs) {
  # to add continuity correction modify the mean as: mean=((n*(m+n+1)-1)/2)
  p.value = stats::pnorm(q=T.obs, mean=((n*(m+n+1))/2), sd = sqrt(n*m*(n+m+1)/12), lower.tail = F)
  return(p.value)
}



# Returns vector of values of the LMPI T3 statistic for each test observation
stat.T3 <- function(Z, m) {
  N = length(Z)
  n = N-m
  R = rank(Z)[(m+1):N]-1
  R2 = R^2
  return(R2+R)
}

# Returns the approximated (1-alpha)-quantile of LMPI T3 based on asymptotic normal approximation
asymptotic.critical.T3 <- function(m, n, alpha) {
  N = m+n
  lambda = m/N
  theta = (4*(2-lambda))/(3*(1-lambda)*lambda)
  variance = 64/(45*N*lambda^3*(1-lambda)^3)
  critical.value = stats::qnorm(alpha, mean=(choose(n,2)*choose(m,2))/N*theta+n*(2*n^2-3*n+1)/6+n*(n-1)/2,
                                sd = (choose(n,2)*choose(m,2))/N*sqrt(variance), lower.tail = F)
  return(critical.value)
}

# Returns the approximated pvalue of theLMPI T3 statistic based on asymptotic normal approximation
asymptotic.pvalue.T3 <- function(m, n, T.obs) {
  N = m+n
  lambda = m/N
  theta = (4*(2-lambda))/(3*(1-lambda)*lambda)
  variance = 64/(45*N*lambda^3*(1-lambda)^3)
  p.value = stats::pnorm(q=T.obs, mean=(choose(n,2)*choose(m,2))/N*theta+n*(2*n^2-3*n+1)/6+n*(n-1)/2,
                         sd = (choose(n,2)*choose(m,2))/N*sqrt(variance), lower.tail = F)
  return(p.value)
}



# Returns vector of values of the Fisher statistic for each test observation
stat.Fisher <- function(Z, m) {
    N = length(Z)
    n = N-m
    S_X = Z[1:m]
    S_Y = Z[(m+1):N]
    pval = sapply(1:n, function(i) (1+sum(S_X >= S_Y[i]))/(m+1))
    R = -2 * log(pval)
    return(R)
}

# Returns the approximated (1-alpha)-quantile of the adjusted Fisher statistic based on asymptotic chi-squared approximation
asymptotic.critical.Fisher <- function(m, n, alpha) {
    gamma = n/m
    critical.value = sqrt(1+gamma) * stats::qchisq(p=1-alpha, df=2*n) - 2 * (sqrt(1+gamma)-1) * n
    return(critical.value)
}

# Returns the approximated pvalue of the adjusted Fisher statistic based on asymptotic chi-squared approximation
asymptotic.pvalue.Fisher <- function(m, n, T.obs) {
  gamma = n/m
  p.value = sqrt(1+gamma) * stats::pchisq(q=T.obs, df=2*n, lower.tail = F) - 2 * (sqrt(1+gamma)-1) * n
  return(p.value)
}


# Returns the (1-alpha)-quantile of the stat.func statistic obtained via permutation
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

# Returns the p-value for the global null obtained via permutation
compute.perm.pval <- function(T.obs, m, n, stat.func, B=10^3, seed=123) {
    set.seed(seed)

    T.v = foreach::foreach(b = 1:B, .combine=cbind) %dopar% {
        N=m+n
        Z = stats::runif(N)
        T = sum(stat.func(Z, m))
        return(T)
    }

    # Compute the permutation p-value
    pval = (1+sum(T.v <= T.obs)) / (1 + length(T.v))

    return(pval)
}

compute.global.pvalue <- function(T.obs, m, n, stat.func, asymptotic.pvalue.func, n_perm, B=100, seed) {

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
#' @param m : calibration size
#' @param n : test size
#' @param alpha : significance level
#' @param stat.func : test statistic of which compute critical values
#' @param asymptotic.critical.func : asymptotic approximation of \code{stat.func}
#' @param n_perm : if \eqn{min(m,n)\leq n_perm} critical values will be computed via permutation. Default value is 10
#' @param B : number of permutation to compute critical values. Default value is 10^3
#' @param critical_values : if not \code{NULL}, a vector of precomputed critical values obtained using
#' the permutation distribution of the test statistic
#' @param seed : seed to ensure reproducible results
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
#' @param S_Y : test score vector
#' @param S_X : calibration score vector
#' @param statistic : parameter indicating the local test to be used in closed testing procedure.
#' It can be either \eqn{T_2, T_3} or adjusted Fisher test
#' @param alpha : significance level
#' @param n_perm : if \eqn{min(m,n)\leq n_perm} critical values will be computed via permutation. Default value is 10
#' @param B : number of permutation to compute critical values. Default value is 10^3
#' @param critical_values : if not \code{NULL}, a vector of precomputed critical values obtained using
#' the permutation distribution of the test statistic
#' @param seed : seed to ensure reproducible results
#'
#' @return An integer which is the \eqn{(1 − \alpha)}-confidence lower bound for
#' the number of true discoveries in closed testing procedure using the chosen local test.
#' No selection in the index test set is performed and the lower bound is computed
#' considering all the observations in the test set.
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
    ##if( (pval.global < alpha) && (d>0) ){
    ##    cat(sprintf("STRANGE. pval.global=%.3f, d=%d.\n", pval.global, d))
    ##}
    out = list("lower.bound" = d, "global.p.value" = pval.global)

    return(out)
}
