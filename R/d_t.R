stat.T2 <- function(Z, m) {
    N = length(Z)
    n = N-m
    R = rank(Z)[(m+1):N]-1
    return(R)
}

asymptotic.critical.T2 <- function(m, n, alpha) {
    critical.value = stats::qnorm(alpha, mean=((n*(m+n+1))/2), sd = sqrt(n*m*(n+m+1)/12), lower.tail = F)
    return(critical.value)
}

stat.T3 <- function(Z, m) {
    N = length(Z)
    n = N-m
    R = rank(Z)[(m+1):N]-1
    R2 = R^2
    return(R2+R)
}

asymptotic.critical.T3 <- function(m, n, alpha) {
    N = m+n
    lambda = m/N
    theta = (4*(2-lambda))/(3*(1-lambda)*lambda)
    variance = 64/(45*N*lambda^3*(1-lambda)^3)
    critical.value = stats::qnorm(alpha, mean=(choose(n,2)*choose(m,2))/N*theta+n*(2*n^2-3*n+1)/6+n*(n-1)/2,
                                  sd = (choose(n,2)*choose(m,2))/N*sqrt(variance), lower.tail = F)
    return(critical.value)
}


perm.crit.T <- function(m, n, stat.func, alpha=0.1, B=10^3, seed=123){
    set.seed(seed)

    T.v = foreach::foreach(b = 1:B, .combine=cbind) %dopar% {
        N=m+n
        Z = runif(N)
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


compute.critical.values <- function(m, n, alpha, stat.func, asymptotic.critical.func, n_perm=10, B=10^3, critical_values=NULL, seed=123){

    crit = sapply(1:n, function(h) {
        if(min(m,h)<=n_perm) {
            found.value = FALSE
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
        } else {
            critical.value = asymptotic.critical.func(m, h, alpha)
        }
        return(critical.value)
    })

    return(crit)
}


d_t <- function(S_Y, S_X, statistic="T2", alpha=0.1, n_perm=10, B=10^3, critical_values=NULL, seed=123){

    stopifnot(statistic %in% c("T2", "T3"))

    if(statistic=="T2") {
        stat.func = stat.T2
        asymptotic.critical.func = asymptotic.critical.T2
    } else if (statistic=="T3") {
        stat.func = stat.T3
        asymptotic.critical.func = asymptotic.critical.T3
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

    return(d)
}
