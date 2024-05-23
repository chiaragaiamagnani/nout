


stat.G.j = function(j, g, N, m, B){

  n = N-m
  stopifnot("j must be less than or equal to N" = j<=N)
  stopifnot("j must be positive" = j>0)

  g.vector = vector(length=B)

  for (b in 1:B) {
    U <- sort(stats::runif(N), decreasing = FALSE)
    U.test <- U[j]

    g.vector[b] <- sapply(X = U.test, FUN = g)
  }

  bar.g.j <- mean(g.vector)

  return(bar.g.j)
}




calc.stat.G = function(B, Z, m, g){

  N = length(Z)
  n = N-m
  R = rank(Z)[(m+1):N] #-1 ?

  if(n==1){
    T.G_j = stat.G.j(R)
  } else {
    T.G_j <- mapply(R, FUN = function(j) stat.G.j(j, g = g, N = N, m = m, B = B))
  }

  T_hat.g <- sum(T.G_j)

  return(T_hat.g)
}




mean.g = function(N, m, g, B){

  n=N-m

  T.G_j <- as.vector(mapply(1:N, FUN = function(j) stat.G.j(j, g = g, N = N, m = m, B = B)))

  out <- n*mean(T.G_j)

  return(out)
}



var.g = function(N, m, g, B){

  n = N-m

  T.G_j <- as.vector(mapply(1:N, FUN = function(j) stat.G.j(j, g = g, N = N, m = m, B = B)))
  mean_T.G_j <- mean(T.G_j)

  out = (m*n*sum((T.G_j - mean_T.G_j)^2))/(N*(N-1))

  return(out)

}


# Analytical counterpart when g = Lehmann's alternative

beta_moment.k = function(h,N,k){

  stopifnot(k>=2)

  r.seq = seq(from=0,to=k-2)
  num = h+r.seq
  den = N+1+r.seq
  out = prod(num/den)

  return(out)
}



mean.analytical.Lk = function(N,k,n){

  stopifnot(k>=2)

  aN_j = k*mapply(1:N, FUN=function(h) beta_moment.k(h,N=N,k=k))
  out = n*mean(aN_j)
  #out = n*bar_aN.Lehmann.k(N=N,k=k)

  return(out)

}




var.analytical.Lk = function(N,k,n){

  stopifnot(k>=2)

  m=N-n

  aN_j = k*mapply(1:N, FUN=function(h) beta_moment.k(h,N=N,k=k))
  bar.aN = mean(aN_j)
  out = n*m*sum((aN_j-bar.aN)^2)/(N*(N-1))

  return(out)

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



# Density function of some of the Lehmann's alternatives
## WMW k=2 (classic WMW)
g2 = function(x, k=2){
  out = ifelse(x<1 & x>0, k*x^(k-1), 0)
  return(out)
}

rg2 = function(n, rg_null, k=2){
  out = replicate(n, max(rg_null(x=k)))
  return(out)
}

## WMW k=3
g3 = function(x, k=3){
  out = ifelse(x<1 & x>0, k*x^(k-1), 0)
  return(out)
}

rg3 = function(n, rg_null, k=3){
  out = replicate(n, max(rg_null(x=k)))
  return(out)
}


## WMW k=10
g10 = function(x, k=10){
  out = ifelse(x<1 & x>0, k*x^(k-1), 0)
  return(out)
}

rg10 = function(n, rg_null, k=10){
  out = replicate(n, max(rg_null(x=k)))
  return(out)
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


asymptotic.critical.Fisher <- function(m, n, alpha=0.1) {
  gamma = n/m
  critical.value = sqrt(1+gamma) * stats::qchisq(p=1-alpha, df=2*n) - 2 * (sqrt(1+gamma)-1) * n
  return(critical.value)
}


asymptotic.critical.Tk <- function(m, n, k, alpha=0.1) {

  stopifnot(k>=1)

  mean.Tk = compute_mean_exact_Tk(m=m,n=n,k=k)
  variance.Tk = n*m*(m+n)/12

  critical.value = stats::qnorm(alpha, mean=mean.Tk, sd = sqrt(variance.Tk), lower.tail = F)

  return(critical.value)
}


compute.critical.values <- function(m, n, alpha, k=NULL, n_perm=10, B=10^3, critical_values=NULL, seed=123){

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

        critical.value = perm.crit.T(m=m, n=h, k=k, alpha=alpha, B=B, seed=seed)
      }
    }
    # For large values of m or n compute critical values using the asymptotic approximation of the test statistic
    else {
      if(is.null(k)){
        critical.value = asymptotic.critical.Fisher(m=m, n=h, alpha=alpha)
      } else {
        critical.value = asymptotic.critical.Tk(m=m, n=h, k=k, alpha=alpha)
      }
    }
    return(critical.value)
  })

  return(crit)
}


