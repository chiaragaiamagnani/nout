
# Utilities that are used to compute variance and mean of Tk in the asymptotic normal approximation


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




Eh0h20.theory = function(k,r=0,s=0,lambda){
  c.00 = factorial(k)^4/(k^2*(1-lambda)^(2*k-2))
  res = ((k-1)^2+2*(k-1))*E2(k=k)+E1(k=k)
  out = c.00*res

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


#' compute_theta.Tk.tilde
#'
#' @param m : calibration sample size
#' @param n : test sample size
#'  @param k : order of the LMPI test statistic
#' @param lambda : proportion of calibration observations in the pooled sample of calibration and test set
#'
#' @return It returns the expectation of LMPI \eqn{T_k} in the limit,
#' computed according to the Central Limit Theorem for \eqn{U}-statistics.
#'
compute_theta.Tk.tilde = function(m,n,k,lambda){
  coef = factorial(k)^2/(k+1)
  somma = sum(sapply(0:(k-1), function(r) choose(k,r)/(lambda^r*(1-lambda)^(k-r-1))))
  out = coef*somma
  return(out)
}



#' compute_theta.Tk
#'
#' @param m : calibration sample size
#' @param n : test sample size
#' @param k : order of the LMPI test statistic
#'
#' @return It returns the expectation of LMPI \eqn{\left[\frac{T_k} -\sum_{i=1}^n(i-1)^k-1\right]{\binom{m,k-1}\binom{n,k-1}}} in the limit,
#' computed according to the Central Limit Theorem for \eqn{U}-statistics.
#'
compute_theta.Tk = function(m,n,k){
  N = n+m
  lambda = m/N

  if(k==1){
    theta = (n*(m+n+1))/2
    #theta = 0.5*(choose(m,k)*choose(n,k)/N) + (n*(n-1))/2

  } else if (k==2) {

    theta = ((4*(2-lambda))/(3*(1-lambda)*lambda))*(choose(m,k)*choose(n,k)/N) + n*(2*n^2-3*n+1)/6 + n*(n-1)/2

  } else {

    coef = factorial(k)^2/(k+1)
    somma = sum(sapply(0:(k-1), function(r) choose(k,r)/(lambda^r*(1-lambda)^(k-r-1))))

    # E[Tk] = E[Tk_tilde]*(choose(m,k)*choose(n,k)/N) + C, dove C = sum(((1:n)-1)^k)
    theta = (coef*somma)*(choose(m,k)*choose(n,k)/N^(k-1)) + sum(((1:n)-1)^k)
  }

  return(theta)
}



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

    # variance = (m*n*(N + 1) / (m^2*n^2*12))*(choose(m,k)*choose(n,k)/N)^2
    variance = n*m*(n+m+1)/12

  } else if (k==2) {

    variance = (64/(45*N*(1-lambda)^3*lambda^3))*(choose(m,k)*choose(n,k)/N)^2

  } else {

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

    theta2 = (compute_theta.Tk.tilde(m=m,n=n,k=k,lambda=lambda))^2
    z10 = Ehh1-theta2
    z01 = Ehh2-theta2
    variance = (k^2/N*(z10/lambda+z01/(1-lambda)))*(choose(m,k)*choose(n,k)/N^(k-1))^2

  }

  return(variance)

}

