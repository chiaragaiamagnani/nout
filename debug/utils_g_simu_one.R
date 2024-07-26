


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

sim_pow = function(B, B_MC, m, n, theta, rg_null, rg, g, g_hat, g_hat_monotone, alpha, q=0.05){

  m = as.double(m)
  n = as.double(n)
  n.out = as.double(round(theta*n))
  n.in = n-n.out
  N = as.double(m+n)

  stats_G_N_oracle =  apply(replicate(B_MC, sapply(X=sort(stats::runif(N)), FUN=g)), 1, mean)
  stats_G_N_KDE = apply(replicate(B_MC, sapply(X=sort(stats::runif(N)), FUN=g_hat)), 1, mean)
  stats_G_N_mono = apply(replicate(B_MC, sapply(X=sort(stats::runif(N)), FUN=g_hat_monotone)), 1, mean)

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

    # Compute test statistics
    T_oracle <-  as.double(calc.stat.G(Z=Z, m=m, stats_G_vector=stats_G_N_oracle))
    T_G_KDE <-  as.double(calc.stat.G(Z=Z, m=m, stats_G_vector=stats_G_N_KDE))
    T_G_mono <-  as.double(calc.stat.G(Z=Z, m=m, stats_G_vector=stats_G_N_mono))
    T_WMW <-  as.double(sum(stat.Tk(Z=Z, m=m, k=1)))
    T_Fisher <-  as.double(sum(stat.Fisher(Z=Z, m=m)))

    res = matrix(data=cbind("T_oracle"=T_oracle, "T_WMW"=T_WMW,
                            "T_Fisher"=T_Fisher, "T_G_KDE"=T_G_KDE, "T_G_mono"=T_G_mono
    ), ncol=5, byrow = T)
    return(res)
  })

  rownames(stats) = c("T_oracle", "T_WMW", "T_Fisher","T_G_KDE", "T_G_mono")

  T_oracle = stats["T_oracle",]
  T_G_KDE = stats["T_G_KDE",]
  T_G_mono= stats["T_G_mono",]
  T_WMW = stats["T_WMW",]
  T_Fisher = stats["T_Fisher",]

  # Compute critical values
  crit_oracle <- as.double(stats::qnorm(alpha, mean=meanG(n=n, stats_G_vector=stats_G_N_oracle),
                              sd = sqrt(varG(n=n, m=m, stats_G_vector=stats_G_N_oracle)), lower.tail = F))
  crit_KDE <-  as.double(stats::qnorm(alpha, mean=meanG(n=n, stats_G_vector=stats_G_N_KDE),
                            sd = sqrt(varG(n=n, m=m, stats_G_vector=stats_G_N_KDE)), lower.tail = F))
  crit_mono <-  as.double(stats::qnorm(alpha, mean=meanG(n=n, stats_G_vector=stats_G_N_mono),
                                 sd = sqrt(varG(n=n, m=m, stats_G_vector=stats_G_N_mono)), lower.tail = F))
  crit_WMW <-  as.double(compute.critical.value.global(m=m, n=n, alpha=alpha, local.test="wmw", k=1, n_perm=10, B=10^3, seed=123))
  crit_Fisher <-  as.double(compute.critical.value.global(m=m, n=n, alpha=alpha, local.test="fisher", k=NULL, n_perm=10, B=10^3, seed=123))


  # Perform the tests
  rej.oracle = ifelse(T_oracle >= crit_oracle, 1, 0)
  rej.KDE = ifelse(T_G_KDE >= crit_KDE, 1, 0)
  rej.mono = ifelse(T_G_mono >= crit_mono, 1, 0)
  rej.WMW = ifelse(T_WMW >= crit_WMW, 1, 0)
  rej.Fisher = ifelse(T_Fisher >= crit_Fisher, 1, 0)

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
               std.error_Fisher = std_err.Fisher
  )

  return(out)

}











sim_d = function(B, B_MC, cons=TRUE, m, n, theta, rg_null, rg, g, g.hat, g_hat_mono, mono, alpha, q=0.05){

  m = as.double(m)
  n = as.double(n)
  n.out = as.double(round(theta*n))
  n.in = n-n.out
  N = as.double(m+n)

  tic()
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

    if(is.null(mono)){
      decr=FALSE
    } else{
      decr=ifelse(mono=="increasing", FALSE, TRUE)
    }

    d_WMW = find_d(X=X, Y=Y, local.test = "wmw", k=1, alpha=alpha, B=B_MC, B_MC=B_MC)$lower.bound
    d_oracle = find_d(X=X, Y=Y, local.test = "g", monotonicity=mono, g.hat=g, alpha=alpha, B=B_MC, B_MC=B_MC)$lower.bound
    if(cons){
      d_estG <- d_G_cons2(S_X=X, S_Y=Y, g.hat=g.hat, alpha=alpha, n_perm=0, B=B_MC, B_MC=B_MC)$lower.bound # è pensato per usare KDE per g e shortcut esatto
    } else{
      d_estG <- d_G_monotone2(S_X=X, S_Y=Y, g.hat=g.hat, decr=decr, alpha=alpha, n_perm=0, B=B_MC, B_MC=B_MC)$lower.bound # è pensato per usare KDE per g e shortcut esatto
    }
    d_estGmono <- d_G_monotone2(S_X=X, S_Y=Y, g.hat=g_hat_mono, decr=decr, alpha=alpha, n_perm=0, B=B_MC, B_MC=B_MC)$lower.bound # è pensato per usare monotone estimate per g e shortcut esatto
    d_fisher <- find_d(X=X, Y=Y, local.test = "fisher", alpha=alpha, B=B_MC, B_MC=B_MC)$lower.bound

    res = matrix(data=cbind("d_WMW"=d_WMW,
                            "d_fisher"=d_fisher,
                            "d_oracle"=d_oracle,
                            "d_estG"=d_estG,
                            "d_estGmono"=d_estGmono
    ), ncol=5, byrow = T)
    return(res)
  })
  toc()

  rownames(stats) = c("d_WMW", "d_fisher", "d_oracle", "d_estG", "d_estGmono")

  d_WMW = stats["d_WMW",]
  d_fisher = stats["d_fisher",]
  d_oracle = stats["d_oracle",]
  d_estG = stats["d_estG",]
  d_estGmono = stats["d_estGmono",]

  power.WMW = mean(d_WMW>0)
  std_err.WMW = (sd(d_WMW>0)/sqrt(B))*qnorm(1-q)

  d_mean.WMW = mean(d_WMW)
  d_std_err.WMW = (sd(d_WMW)/sqrt(B))*qnorm(1-q)

  power.fisher = mean(d_fisher>0)
  std_err.fisher = (sd(d_fisher>0)/sqrt(B))*qnorm(1-q)

  d_mean.fisher = mean(d_fisher)
  d_std_err.fisher = (sd(d_fisher)/sqrt(B))*qnorm(1-q)

  power.oracle= mean(d_oracle>0)
  std_err.oracle = (sd(d_oracle>0)/sqrt(B))*qnorm(1-q)

  d_mean.oracle = mean(d_oracle)
  d_std_err.oracle = (sd(d_oracle)/sqrt(B))*qnorm(1-q)

  power.estG = mean(d_estG>0)
  std_err.estG = (sd(d_estG>0)/sqrt(B))*qnorm(1-q)

  d_mean.estG = mean(d_estG)
  d_std_err.estG = (sd(d_estG)/sqrt(B))*qnorm(1-q)

  power.estGmono = mean(d_estGmono>0)
  std_err.estGmono = (sd(d_estGmono>0)/sqrt(B))*qnorm(1-q)

  d_mean.estGmono = mean(d_estGmono)
  d_std_err.estGmono = (sd(d_estGmono)/sqrt(B))*qnorm(1-q)

  out = tibble(n = n,
               power_oracle = power.oracle,
               std.error_oracle = std_err.oracle,
               d_oracle = d_mean.oracle,
               d.std.error_oracle = d_std_err.oracle,
               power_estG = power.estG,
               std.error_estG = std_err.estG,
               d_estG = d_mean.estG,
               d.std.error_estG = d_std_err.estG,
               power_estGmono = power.estGmono,
               std.error_estGmono = std_err.estGmono,
               d_estGmono = d_mean.estGmono,
               d.std.error_estGmono = d_std_err.estGmono,
               power_WMW = power.WMW,
               std.error_WMW = std_err.WMW,
               d_WMW = d_mean.WMW,
               d.std.error_WMW = d_std_err.WMW,
               power_Fisher = power.fisher,
               std.error_Fisher = std_err.fisher,
               d_Fisher = d_mean.fisher,
               d.std.error_Fisher = d_std_err.fisher
  )

  return(out)

}

