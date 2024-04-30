
# --------------------------------------------------------- MEAN T_k TILDE ---------------------------------------------------------


# This function computes the approximated mean of Tk tilde applying Theorem 1
compute_approx_mean.Tk.tilde = function(m,n,k){
  N = n+m
  lambda = m/N

  if(k==1){ # WMW
    theta = 0.5

  } else if (k==2) {

    theta = (4*(2-lambda))/(3*(1-lambda)*lambda)

  } else if(k>2) {

    coef = factorial(k)^2/(k+1)
    somma = sum(sapply(0:(k-1), function(r) choose(k,r)/(lambda^r*(1-lambda)^(k-r-1))))
    theta = (coef*somma)
  }

  return(theta)
}



# It computes the exact mean of Tk tilde using the inverse transformation
compute_mean_exact_Tk.tilde = function(m,n,k){

  N=n+m

  Ck = calc.const.Tk(n=n,k=k)
  mean = compute_mean_exact_Tk(m=m,n=n,k=k)

  if(k<4){
    out = ((mean-Ck)*N^(k-1))/(choose(n,k)*choose(m,k))
  } else {
    out = NA
  }

  return(out)
}









# --------------------------------------------------------- MEAN T_k ---------------------------------------------------------



From_Tk.tilde_To_Tk_mean = function(mean.tilde,m,n,k){

  stopifnot(k>=1)

  N=n+m

  mean = mean.tilde*(choose(m,k)*choose(n,k)/N^(k-1))

  return(mean)
}



# This function computes the approximated mean of Tk applying Theorem 1
compute_approx_mean.Tk = function(m,n,k){

  N = n+m

  # E[Tk] = E[Tk_tilde]*(choose(m,k)*choose(n,k)/N) + C, dove C = sum(((1:n)-1)^k)
  Ck = calc.const.Tk(n=n,k=k)
  approx_mean.Tk.tilde=compute_approx_mean.Tk.tilde(m=m,n=n,k=k)
  mu = From_Tk.tilde_To_Tk_mean(mean.tilde=approx_mean.Tk.tilde,m=m,n=n,k=k) + Ck # sum(((1:n)-1)^k)

  return(mu)
}



# Compute mu_k (= expected value of R_1^k) according to Lemma 2.
# This algorithm has been formally proved for k=1,2,3 but not for further values of k.
compute.mu_k = function(m,n,k){

  l_v = 0:k
  c.binom_v = choose(k,l_v)

  E = vector()

  if(k==1){

    k.bar_n = min(n,k)
    k.bar_m = min(m,k)

    E0 = sum(sapply(1:k.bar_m, function(dj){ prod((m-dj+1):m) / (dj+1)}))
    Ek = sum(sapply(1:k.bar_n, function(di){ prod(((n-di):(n-1))) / (di+1)}))
    E = c(E0,Ek)

  } else if(k==2){
    k.bar_n = min(n,k)
    k.bar_m = min(m,k)

    dj_v = 1:k.bar_m
    di_v = 1:k.bar_n

    c.dj_v = c(1,choose(k.bar_m, dj_v[2:k.bar_m]))
    c.di_v = c(1,choose(k.bar_n, di_v[2:k.bar_n]))

    E0 = sum(sapply(1:k.bar_m, function(dj){c.dj_v[dj] * prod(((m-dj+1):m)) / (dj+1)}))
    Ek = sum(sapply(1:k.bar_n, function(di){c.di_v[di] * prod(((n-di):(n-1))) / (di+1)}))

    for(l in 1:(k-1)){
      k.bar_n = min(n,l)
      k.bar_m = min(m,k-l)

      dj_v = 1:(k.bar_m) # k.bar_n
      di_v = 1:k.bar_n

      if((k.bar_m-l)==1){
        c.dj_v = choose(k.bar_m, dj_v)
      } else {
        c.dj_v = c(1,choose(k.bar_m, dj_v[2:(k.bar_m)]))
      }

      if(k.bar_n==1){
        c.di_v = choose(k.bar_n, di_v)
      } else {
        c.di_v = c(1,choose(k.bar_n, di_v[2:k.bar_n]))
      }

      E[l] = sum(sapply(dj_v, function(dj){
        sum(sapply(di_v, function(di){
          c.dj_v[dj]*c.di_v[di]*prod(((m-dj+1):m))*prod(((n-di):(n-1))) / (di+dj+1)
        }))
      }))
    }

    E = c(E0,E,Ek)

  } else {

    k.bar_n = min(n,k)
    k.bar_m = min(m,k)

    dj_v = 1:k.bar_m
    di_v = 1:k.bar_n

    c.dj_v = c(1,choose(k.bar_m, dj_v[2:k.bar_m]))
    c.di_v = c(1,choose(k.bar_n, di_v[2:k.bar_n]))

    E0 = sum(sapply(1:k.bar_m, function(dj){c.dj_v[dj] * prod(((m-dj+1):m)) / (dj+1)}))
    Ek = sum(sapply(1:k.bar_n, function(di){c.di_v[di] * prod(((n-di):(n-1))) / (di+1)}))


    for(l in 1:(k-1)){

      k.bar_n = min(n,l)
      k.bar_m = min(m,k-l)

      dj_v = 1:k.bar_m
      di_v = 1:k.bar_n

      if(k.bar_m==1){
        c.dj_v = choose(k.bar_m, dj_v)
      } else {
        c.dj_v = c(1,choose(k.bar_m, dj_v[2:k.bar_m]))
      }

      if(k.bar_n==1){
        c.di_v = choose(k.bar_n, di_v)
      } else {
        c.di_v = c(1,choose(k.bar_n, di_v[2:k.bar_n]))
      }

      E[l] = sum(sapply(dj_v, function(dj){
        sum(sapply(di_v, function(di){
          c.dj_v[dj]*c.di_v[di]*prod(((m-dj+1):m))*prod(((n-di):(n-1))) / (di+dj+1)
        }))
      }))

    }

    E = c(E0,E,Ek)
  }

  mu_k = sum(c.binom_v*E)

  return(mu_k)
}


# Compute the exact mean of Tk according to Lemma 1
compute_mean_exact_Tk = function(m,n,k){

  if(k==1){
    out = n*(compute.mu_k(m=m,n=n,k=1))
  } else if(k==2){
    out = n*(compute.mu_k(m=m,n=n,k=1) + compute.mu_k(m=m,n=n,k=2))
  } else if(k==3){
    out = n*(2*compute.mu_k(m=m,n=n,k=1) + 3*compute.mu_k(m=m,n=n,k=2) + compute.mu_k(m=m,n=n,k=3))
  } else {
    out = NA
  }

  return(out)
}




# --------------------------------------------------------- VARIANCE T_k TILDE ---------------------------------------------------------

E3 = function(k){
  a1 = numeric()
  a2 = numeric()

  for(h in 0:(k-1)){

    aa1 = numeric()

    for(t in 0:(2*k-2*h)){
      aa1[t+1] = choose(2*k-2*h, t)*(-1)^t/(2*h+t+1)
    }
    a1[h+1] = (factorial(k-1))^2/(factorial(k-h)*factorial(h))^2 * sum(aa1)
  }
  A1 = sum(a1)


  for(h in 0:(k-2)){

    aa2 = numeric()

    for(s in (h+1):(k-1)) {

      c2 = factorial(k-1)^2/(factorial(k-h)*factorial(h)*factorial(k-s)*factorial(s))

      aaa2 = numeric()
      for(t in 0:(2*k-h-s)){
        aaa2[t+1] = choose(2*k-h-s, t)*(-1)^t/(h+s+t+1)
      }
      aa2[s-h]=2*c2*sum(aaa2)
    }
    a2[h+1]=sum(aa2)
  }
  A2 = sum(a2)

  out = A1+A2

  return(out)
}



E2 = function(k){
  out = 1/((k+1)^2)
  return(out)
}



E1 = function(k){
  out = 1/(2*k+1)
  return(out)
}


E0 = function(k){
  a1 = numeric()

  for(h in 0:(k-1)){

    aa1 = numeric()

    for(t in 0:(k-h)){
      aa1[t+1] = choose(k-h, t)*(-1)^t/(k+h+t+1)
    }
    a1[h+1] = sum(aa1)/(factorial(k-h)*factorial(h))
  }
  A1 = sum(a1)

  return(factorial(k-1)*A1)
}


Ehrhr.theory = function(k,r,lambda){

  if(r==0){
    c.out = factorial(k)^4/(k^2*(1-lambda)^(2*k-2))
    out = c.out*k^2*E3(k=k)
  }
  else{
    c.out = (choose(k,r)*factorial(k-r)*factorial(k-r-1)*factorial(r)^2 / (lambda^r*(1-lambda)^(k-r-1)))^2
    res = (choose(k-1,k-r-1)*k*choose(k-1,r))^2*E3(k=k)+
      2*choose(k-1,k-r-1)*k*choose(k-1,r)*choose(k-1,k-r)*k*choose(k-1,r)*E2(k=k)+
      (choose(k-1,k-r)*k*choose(k-1,r))^2*E2(k=k)
    out = c.out*res
  }

  return(out)
}



Eh0hr.theory = function(k,r,lambda){

  if(r==0){
    c.out = factorial(k)^4/(k^2*(1-lambda)^(2*k-2))
    out = c.out*k^2*E3(k=k)
  }
  else{
    c.out = (choose(k,r)*factorial(k-r)*factorial(k-r-1)*factorial(r)^2 / (lambda^r*(1-lambda)^(k-r-1))) *
      factorial(k)^2/(k*(1-lambda)^(k-1))
    res = (choose(k-1,k-r-1)*k^2*choose(k-1,r))*E3(k=k)+
      (choose(k-1,k-r)*k^2*choose(k-1,r))*E2(k=k)
    out = c.out*res
  }

  return(out)
}



Ehrhs.theory = function(k,r,s,lambda){

  c.r = (choose(k,r)*factorial(k-r)*factorial(k-r-1)*factorial(r)^2 / (lambda^r*(1-lambda)^(k-r-1)))
  c.s = (choose(k,s)*factorial(k-s)*factorial(k-s-1)*factorial(s)^2 / (lambda^s*(1-lambda)^(k-s-1)))

  res = (choose(k-1,k-r-1)*k^2*choose(k-1,r)*choose(k-1,k-s-1)*choose(k-1,s))*E3(k=k)+
    (choose(k-1,k-r-1)*k^2*choose(k-1,r)*choose(k-1,k-s)*choose(k-1,s))*E2(k=k)+
    (choose(k-1,k-r)*k^2*choose(k-1,r)*choose(k-1,k-s)*choose(k-1,s))*E2(k=k)+
    (choose(k-1,k-r)*k^2*choose(k-1,r)*choose(k-1,k-s-1)*choose(k-1,s))*E2(k=k)

  out = c.r*c.s*res

  return(out)
}


Eh0h2r.theory = function(k,r=0,s=NULL,lambda){
  c.0 = factorial(k)^2/(k*(1-lambda)^(k-1))
  c.r = choose(k,r)*factorial(k-r)*factorial(k-r-1)*factorial(r)^2/(lambda^r*(1-lambda)^(k-r-1))
  e2 = E2(k=k)
  res = (k-1)*choose(k,k-r)*(k*choose(k-1,r)-choose(k-1,r)-(k-1)*choose(k-2,r-1))*e2 +
    (k-1)*choose(k,k-r)*choose(k-1,r)*e2 +
    (k-1)^2*choose(k,k-r)*choose(k-2,r-1)*e2 +
    choose(k,k-r)*(k*choose(k-1,r)-choose(k-1,r)-(k-1)*choose(k-2,r-1))*e2 +
    choose(k,k-r)*choose(k-1,r)*E1(k=k) +
    choose(k,k-r)*(k-1)*choose(k-2,r-1)*E0(k=k)

  out = c.0*c.r*res

  return(out)
}

Eh0h20.theory = function(k,r=0,s=0,lambda){
  c.00 = factorial(k)^4/(k^2*(1-lambda)^(2*k-2))
  res = ((k-1)^2+2*(k-1))*E2(k=k)+E1(k=k)
  out = c.00*res

  return(out)
}

Ehrh2r.theory = function(k,r,s=NULL,lambda){

  c.r = (choose(k,r)*factorial(k-r)*factorial(k-r-1)*factorial(r)^2 / (lambda^r*(1-lambda)^(k-r-1)))

  e2 = E2(k=k)
  res = (choose(k,k-r)*choose(k-1,r))^2 * E1(k=k)+
    2*(choose(k,k-r)^2*(k-1)*choose(k-2,r-1)*choose(k-1,r)) * E0(k=k)+
    (choose(k,k-r)*(k-1)*choose(k-2,r-1))^2 * E3(k=k)+
    (choose(k,k-r)^2*(k-1)*((k-1)*(choose(k-1,r)-choose(k-2,r-1))))*e2+
    (choose(k,k-r)^2*(k-1)*choose(k-2,r-1)*((k-1)*(choose(k-1,r)-choose(k-2,r-1))))*e2+
    (choose(k,k-r)^2*(k)*choose(k-1,r)*((k-1)*(choose(k-1,r)-choose(k-2,r-1))))*e2

  out = c.r^2*res

  return(out)
}



Ehrh2s.theory = function(k,r,s,lambda){

  c.r = (choose(k,r)*factorial(k-r)*factorial(k-r-1)*factorial(r)^2 / (lambda^r*(1-lambda)^(k-r-1)))
  c.s = (choose(k,s)*factorial(k-s)*factorial(k-s-1)*factorial(s)^2 / (lambda^s*(1-lambda)^(k-s-1)))

  e2 = E2(k=k)
  e0 = E0(k=k)

  res = (choose(k,k-r)*choose(k-1,r)*choose(k,k-s)*choose(k-1,s)) * E1(k=k) +

    (choose(k,k-r)*choose(k-1,r)*choose(k,k-s)*(k-1)*choose(k-2,s-1)) * e0 +

    (choose(k,k-r)*choose(k-1,r)*choose(k,k-s)*(k-1)*(choose(k-1,s)-choose(k-2,s-1))) * e2 +

    (choose(k,k-r)*(k-1)*choose(k-2,r-1)*choose(k,k-s)*choose(k-1,s)) * e0 +

    (choose(k,k-r)*(k-1)^2*choose(k-2,r-1)*choose(k,k-s)*choose(k-2,s-1)) * E3(k=k) +

    (choose(k,k-r)*(k-1)^2*choose(k-2,r-1)*choose(k,k-s)*(choose(k-1,s)-choose(k-2,s-1))) * e2 +

    (choose(k,k-r)*(k-1)*(choose(k-1,r)-choose(k-2,r-1))*choose(k,k-s)*k*(choose(k-1,s))) * e2

  out = c.r*c.s*res

  return(out)
}



# Compute variance of Tk tilde according to Theorem 1
compute_var.Tk.tilde = function(k, m, n){

  N = n+m
  lambda = m/N

  if(k==1){ # WMW

    # variance = (m*n*N / (m^2*n^2*12))*(choose(m,k)*choose(n,k))^2
    variance = N/(12*m*n)

  } else if (k==2) { # T3

    variance = (64/(45*N*(1-lambda)^3*lambda^3))

  } else if(k>2) {

    Ehh1 = sum(sapply(0:(k-1), function(r) Ehrhr.theory(k=k, r=r, lambda=lambda))) +
      2 * sum(sapply(1:(k-2), function(r){
        sum(sapply((r+1):(k-1), function(s) Ehrhs.theory(k=k, r=r, s=s, lambda=lambda)))
      })) +
      2 * sum(sapply(1:(k-1), function(r) Eh0hr.theory(k=k, r=r, lambda=lambda)
      ))

    Ehh2 = Eh0h20.theory(k=k, r=0, lambda=lambda) +
      sum(sapply(1:(k-1), function(r) Ehrh2r.theory(k=k, r=r, lambda=lambda))) +
      2 * sum(sapply(1:(k-2), function(r){
        sum(sapply((r+1):(k-1), function(s) Ehrh2s.theory(k=k, r=r, s=s, lambda=lambda)))
      })) +
      2 * sum(sapply(1:(k-1), function(r) Eh0h2r.theory(k=k, r=r, lambda=lambda)
      ))

    theta2 = (compute_approx_mean.Tk.tilde(m=m,n=n,k=k))^2
    z10 = Ehh1-theta2
    z01 = Ehh2-theta2
    variance = (k^2/N*(z10/lambda+z01/(1-lambda)))

  }

  return(variance)

}




# --------------------------------------------------------- VARIANCE T_k ---------------------------------------------------------


#' compute_variance.Tk
#'
#' @param m : calibration sample size
#' @param n : test sample size
#' @param k : order of the LMPI test statistic
#'
#' @return It returns the variance of LMPI \eqn{T_k} in the limit,
#' computed according to the Central Limit Theorem for \eqn{U}-statistics.
#'
compute_variance.Tk = function(k, m, n){

  N = n+m
  lambda = m/N

  if(k==1){

    # variance = (m*n*N / (m^2*n^2*12))*(choose(m,k)*choose(n,k))^2
    variance = n*m*N/12

  } else {

    variance = compute_var.Tk.tilde(m=m,n=n,k=k)*(choose(m,k)*choose(n,k)/N^(k-1))^2

  }

  return(variance)

}



From_Tk.tilde_To_Tk_variance = function(variance.tilde,m,n,k){

  stopifnot(k>=1)

  N=n+m

  variance = variance.tilde*(choose(m,k)*choose(n,k)/N^(k-1))^2

  return(variance)
}





# ----------------------------------------------------------- GENERATE DATA T_k -----------------------------------------------------------


gen.data <- function(m,n) {
  Z <- stats::rnorm(m+n)
  return(Z)
}

stat.Tk <- function(Z, m, k) {

  N = length(Z)
  n = N-m
  R = rank(Z)[(m+1):N]-1

  if(n==1){
    Tk_i = sum(sapply(1:k, function(h){R^h}))
  }
  else{
    # Tk_i = apply(sapply(1:k, function(h){R^h}), MARGIN=1, FUN = sum)
    Tk_i = apply(sapply(1:k, function(l){R+l-1}), MARGIN=1, FUN = prod)
  }

  return(Tk_i)
}


# New notation (T1 = WMW)
calc.Tk <- function(Z, m, k) {
  N <- length(Z)
  Tk = sum(stat.Tk(Z=Z,m=m,k=k))
  return(Tk)
}


# New notation (T1 = WMW)
calc.const.Tk = function(n,k){

  if(k==1){ # WMW
    Ck <- n*(n-1)/2
  } else if(k==2){
    Ck <- n*(2*n^2-3*n+1)/6 + n*(n-1)/2
  } else if(k>2){
    Ck <- sum(sapply(1:k, function(h) sum(((1:n)-1)^h)))
  }
  return(Ck)
}


calc.Tk.tilde <- function(Z,m,k) {
  N <- length(Z)
  n <- N-m
  Tk <- calc.Tk(Z=Z,m=m,k=k)

  if(k==1){ # WMW
    Tkt <- Tk - n*(n-1)/2
  } else if(k==2){
    Tkt <- Tk - n*(2*n^2-3*n+1)/6 - n*(n-1)/2
  } else if(k>2){
    Ck = calc.const.Tk(n=n,k=k)
    Tkt <- Tk - Ck
  }

  Tkt <- N^(k-1) * Tkt / (choose(n,k) * choose(m,k))
  return(Tkt)
}





# ----------------------------------------------------- GENERATE DATA sum_{i\in[n]} R_i^k -----------------------------------------------------


# New notation (T1 = WMW)
calc.sumRik <- function(Z, m, k) {
  N <- length(Z)
  sumRk = sum(stat.Rik(Z=Z,m=m,k=k))
  return(sumRk)
}


stat.Rik <- function(Z, m, k) {

  N = length(Z)

  n = N-m
  R = rank(Z)[(m+1):N]-1

  Rk_i=R^k

  return(Rk_i)
}


# New notation (T1 = WMW)
calc.const.sumRk = function(n,k){

  if(k==1){ # WMW
    Ck <- n*(n-1)/2
  } else if(k==2){
    Ck <- n*(2*n^2-3*n+1)/6
  } else if(k==3){
    Ck <- n^2*(n-1)^2/4 - 3*n*(n-1)/2
  } else if(k>3){
    Ck <- sum(((1:n)-1)^k)
  }
  return(Ck)
}


calc.sumRk.tilde <- function(Z,m,k) {
  N <- length(Z)
  n <- N-m
  sumRk <- calc.sumRik(Z=Z,m=m,k=k)

  if(k==1){ # WMW
    sumRkt <- sumRk - n*(n-1)/2
  } else if(k==2){
    sumRkt <- sumRk - n*(2*n^2-3*n+1)/6
  } else if(k>2){
    Ck = calc.const.sumRk(n=n,k=k)
    sumRkt <- sumRk - Ck
  }

  sumRkt <- N^(k-1) * sumRkt / (choose(n,k) * choose(m,k))
  return(sumRkt)
}






