


# New notation (T1 = WMW)
calc.Tk <- function(Z, m, k) {
  m = as.double(m)
  k = as.double(k)
  N <- as.double(length(Z))
  Tk = sum(stat.Tk(Z=Z,m=m,k=k))
  return(Tk)
}





# Compute LMPI test statistic using Shiraishi's method (1985) against Lehmann's alternative

## Test statistic for one single test point
stat.Tk_j_Shiraishi <- function(Z, m, k) {

  stopifnot("k must be greater than or equal to 2" = k>=2)

  m = as.double(m)
  k = as.double(k)
  N = as.double(length(Z))
  n = as.double(N-m)
  R = rank(Z)[(m+1):N]#-1

  if(n==1){
    Tk_i = prod(sapply(0:(k-2), function(l){(R+l)/(N+1+l)}))
  } else {
    Tk_i = apply(sapply(0:(k-2), function(l){(R+l)/(N+1+l)}), MARGIN=1, FUN = prod)
  }
  return(Tk_i)
}


## Test statistic for n test points
calc.Tk_Shiraishi <- function(Z, m, k) {
  stopifnot("k must be greater than or equal to 2" = k>=2)

  m = as.double(m)
  k = as.double(k)
  N <- as.double(length(Z))
  Tk = k*sum(stat.Tk_j_Shiraishi(Z=Z,m=m,k=k))
  return(Tk)
}



# Functions not exported in the R package nout but needed for simulations

stat.Fisher <- function(Z, m) {
  N = length(Z)
  n = N-m
  S_X = Z[1:m]
  S_Y = Z[(m+1):N]
  pval = sapply(1:n, function(i) (1+sum(S_X >= S_Y[i]))/(m+1))
  R = -2 * log(pval)
  return(R)
}


# perm.crit.T <- function(m, n, k=NULL, alpha=0.1, B=10^3, seed=123){
#   set.seed(seed)
#
#   T.v = foreach::foreach(b = 1:B, .combine=cbind) %dopar% {
#     N=m+n
#     Z = stats::runif(N)
#     if(is.null(k)){
#       T = sum(stat.Fisher(Z=Z, m=m))
#     }
#     else {
#       T = sum(stat.Tk(Z=Z, k=k, m=m))
#     }
#   }
#
#   # Empirical quantile with finite-sample adjustment
#   idx = ceiling(B*(1-alpha)*(1+1/B))
#   if(idx <= B) {
#     crit = sort(as.vector(T.v))[idx]
#   } else {
#     crit = Inf
#   }
#
#   return(crit)
# }


asymptotic.critical.Fisher <- function(m, n, alpha=0.1) {
  gamma = n/m
  critical.value = sqrt(1+gamma) * stats::qchisq(p=1-alpha, df=2*n) - 2 * (sqrt(1+gamma)-1) * n
  return(critical.value)
}


asymptotic.critical.Tk <- function(m, n, k, alpha=0.1) {

  stopifnot(k>=1)

  # mean.Tk = compute_mean_exact_Tk(m=m,n=n,k=k)
  mean.Tk = m*n/2 + n*(n-1)/2
  variance.Tk = n*m*(m+n)/12

  critical.value = stats::qnorm(alpha, mean=mean.Tk, sd = sqrt(variance.Tk), lower.tail = F)

  return(critical.value)
}

#
# compute.critical.values <- function(m, n, alpha, k=NULL, n_perm=10, B=10^3, critical_values=NULL, seed=123){
#
#   crit = sapply(1:n, function(h) {
#     # For small values of m and n compute critical values via permutation
#     if(min(m,h)<=n_perm) {
#       found.value = FALSE
#       # In order to avoid repeating computation of precomputed critical values that are saved in "tables" folder
#       if(!is.null(critical_values)) {
#         if(length(critical_values)>=h) {
#           critical.value = critical_values[h]
#           found.value = TRUE
#         }
#       }
#       if(!found.value) {
#         #cat(sprintf("Running permutations...\n"))
#
#         critical.value = perm.crit.T(m=m, n=h, k=k, alpha=alpha, B=B, seed=seed)
#       }
#     }
#     # For large values of m or n compute critical values using the asymptotic approximation of the test statistic
#     else {
#       if(is.null(k)){
#         critical.value = asymptotic.critical.Fisher(m=m, n=h, alpha=alpha)
#       } else {
#         critical.value = asymptotic.critical.Tk(m=m, n=h, k=k, alpha=alpha)
#       }
#     }
#     return(critical.value)
#   })
#
#   return(crit)
# }


compute.critical.value.global <- function(m, n, alpha, local.test, k=NULL, n_perm=10, B=10^3, seed=123){

  # For small values of m and n compute critical values via permutation
  if(min(m,n)<=n_perm) {
    critical.value = perm.crit.T(m=m, n=n, k=k, local.test=local.test, alpha=alpha, B=B, seed=seed)
  }
  # For large values of m or n compute critical values using the asymptotic approximation of the test statistic
  else {
    if(is.null(k)){
      critical.value = asymptotic.critical.Fisher(m=m, n=n, alpha=alpha)
    } else {
      critical.value = asymptotic.critical.Tk(m=m, n=n, k=k, alpha=alpha)
    }
  }

  return(critical.value)
}








# Simulation comparing the power of the oracle test, WMW test, Fisher's method and Shiraishi test

sim_pow_rep = function(B, B_MC, m, n, theta, rg_null, rg, g, constraint="increasing", propF = 0.5, alpha, q=0.05){

  m = as.double(m)
  n = as.double(n)
  n.out = as.double(round(theta*n))
  n.in = n-n.out
  N = as.double(m+n)

  stats <- sapply(1:B, function(b) {

    # Generate scores
    X <- rg_null(m)
    Y.in <- rg_null(n.in)
    if(n.out>0){
      Y.out <- rg(n.out)
      Y <- c(Y.in, Y.out)
    } else {
      Y <- Y.in
    }
    Z <- c(X,Y)
    X1 = sample(X, size = propF*m)
    X2 = setdiff(X, X1)

    g_mono = estimate_g(X1=X1, X2=X2, Y=Y, constraint=constraint, ker="uniform")
    g_KDE = estimate_g(X1=X1, X2=X2, Y=Y, ker="uniform")

    stats_G_N_oracle =  apply(replicate(B_MC, sapply(X=sort(stats::runif(N)), FUN=g)), 1, mean)
    stats_G_N_KDE = apply(replicate(B_MC, sapply(X=sort(stats::runif(N)), FUN=g_KDE)), 1, mean)
    stats_G_N_mono = apply(replicate(B_MC, sapply(X=sort(stats::runif(N)), FUN=g_mono)), 1, mean)

    # Compute test statistics
    T_oracle <- calc.stat.G(Z=Z, m=m, stats_G_vector=stats_G_N_oracle)
    T_G_KDE <- calc.stat.G(Z=Z, m=m, stats_G_vector=stats_G_N_KDE)
    T_G_mono <- calc.stat.G(Z=Z, m=m, stats_G_vector=stats_G_N_mono)
    T_WMW <- calc.Tk(Z=Z, m=m, k=1)
    T_Fisher <- sum(stat.Fisher(Z=Z, m=m))

    crit_oracle <-  as.double(stats::qnorm(alpha, mean=meanG(n=n, stats_G_vector=stats_G_N_oracle),
                                sd = sqrt(varG(n=n, m=m, stats_G_vector=stats_G_N_oracle)), lower.tail = F))
    crit_estG <-  as.double(stats::qnorm(alpha, mean=meanG(n=n, stats_G_vector=stats_G_N_KDE),
                              sd = sqrt(varG(n=n, m=m, stats_G_vector=stats_G_N_KDE)), lower.tail = F))
    crit_estG_mono <-  as.double(stats::qnorm(alpha, mean=meanG(n=n, stats_G_vector=stats_G_N_mono),
                                   sd = sqrt(varG(n=n, m=m, stats_G_vector=stats_G_N_mono)), lower.tail = F))
    crit_WMW <-  as.double(compute.critical.value.global(m=m, n=n, alpha=alpha, local.test="wmw", k=1, n_perm=0, B=B, seed=123))
    crit_Fisher <-  as.double(compute.critical.value.global(m=m, n=n, alpha=alpha, local.test="fisher", k=NULL, n_perm=0, B=B, seed=123))

    rej.oracle = ifelse(T_oracle >= crit_oracle, 1, 0)
    rej.estG = ifelse(T_G_KDE >= crit_estG, 1, 0)
    rej.estG_mono = ifelse(T_G_mono >= crit_estG_mono, 1, 0)
    rej.WMW = ifelse(T_WMW >= crit_WMW, 1, 0)
    rej.Fisher = ifelse(T_Fisher >= crit_Fisher, 1, 0)

    res = matrix(data=cbind("rej_oracle"=rej.oracle, "rej_WMW"=rej.WMW,
                            "rej_Fisher"=rej.Fisher,
                            "rej_G_KDE"=rej.estG ,
                            "rej_G_mono"=rej.estG_mono#, "KDE"=g_KDE, "mono"=g_mono
    ), ncol=5, byrow = T)

    return(res)
  })

  rownames(stats) = c("rej_oracle", "rej_WMW", "rej_Fisher", "rej_G_KDE", "rej_G_mono")


  rej.oracle = unlist(stats[1,])
  rej.WMW = unlist(stats[2,])
  rej.Fisher = unlist(stats[3,])
  rej.KDE = unlist(stats[4,])
  rej.mono = unlist(stats[5,])

  power.oracle = mean(rej.oracle)
  std_err.oracle = (sd(rej.oracle)/sqrt(B))*qnorm(1-q)

  power.KDE = mean(rej.KDE)
  std_err.KDE = (sd(rej.KDE)/sqrt(B))*qnorm(1-q)

  power.mono = mean(rej.mono)
  std_err.mono = (sd(rej.mono)/sqrt(B))*qnorm(1-q)

  power.WMW = mean(rej.WMW)
  std_err.WMW = (sd(rej.WMW)/sqrt(B))*qnorm(1-q)

  power.Fisher = mean(rej.Fisher)
  std_err.Fisher = (sd(rej.Fisher)/sqrt(B))*qnorm(1-q)

  out = tibble(n = n,
               power_oracle = power.oracle,
               std.error_oracle = std_err.oracle,
               power_mono = power.mono,
               std.error_mono = std_err.mono,
               power_KDE = power.KDE,
               std.error_KDE = std_err.KDE,
               power_WMW = power.WMW,
               std.error_WMW = std_err.WMW,
               power_Fisher = power.Fisher,
               std.error_Fisher = std_err.Fisher)


  return(out)

}
