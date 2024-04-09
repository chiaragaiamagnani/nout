
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
#' @param k : order of the LMPI test statistic
#'
#' @return It returns the expectation of LMPI \eqn{T_k} in the limit,
#' computed according to the Central Limit Theorem for \eqn{U}-statistics.
#'
compute_theta.Tk.tilde = function(m,n,k){
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

    # variance = (m*n*(N + 1) / (m^2*n^2*12))*(choose(m,k)*choose(n,k))^2
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

    theta2 = (compute_theta.Tk.tilde(m=m,n=n,k=k))^2
    z10 = Ehh1-theta2
    z01 = Ehh2-theta2
    variance = (k^2/N*(z10/lambda+z01/(1-lambda)))*(choose(m,k)*choose(n,k)/N^(k-1))^2

  }

  return(variance)

}


calc.theory.variance.Tk.tilde = function(k, m, n){

  N = n+m
  lambda = m/N

  if(k==1){ # WMW

    # variance = (m*n*(N + 1) / (m^2*n^2*12))*(choose(m,k)*choose(n,k))^2
    variance = (N+1)/(12*m*n)

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

    theta2 = (compute_theta.Tk.tilde(m=m,n=n,k=k))^2
    z10 = Ehh1-theta2
    z01 = Ehh2-theta2
    variance = (k^2/N*(z10/lambda+z01/(1-lambda)))

  }

  return(variance)

}


# This function computes the variance of the remainder of T2 tilde
calc.theory.variance.remainderT2.tilde = function(m,n){
  N = n+m
  z01_Re = (N^2*(16*m^4 + 16*n^4 + 30*m^3*(-1 + 2*n) +
                   30*m*n^2*(-1 + 2*n) + m^2*(15 - 60*n + 92*n^2)))/(45*m^4*n^4)
  z10_Re = (N^2*(16*m^4 + 16*n^4 + 30*m^3*(-1 + 2*n) +
                   30*m*n^2*(-1 + 2*n) + m^2*(75 + 4*n*(-75 + 83*n))))/(45*m^4*n^4)
  variance.remainder =  4*(z10_Re*N/m+z01_Re*N/n)

  return(variance.remainder)
}

# This function computes the variance of the remainder of T2
# calc.theory.variance.remainderT2 = function(m,n){
#
#   N = n+m
#   variance.remainder =  calc.theory.variance.remainder.tilde(m=m,n=n)/N
#
#   return(variance.remainder)
# }

calc.theory.mean.remainderT3.tilde = function(m,n){

  N = n+m
  lambda = m/N

  mean.remainderT3 = (81*(4 + 3*m))/(lambda^2 * m^2) +
    ( 27*(3 - 2*m + m*n - 2*n))/((1 - lambda) * lambda * m^2 * n) +
    ( 36*(-3 + n)) / ((1 - lambda)^2*m*n) +
    (81*(-4 + m + 2*n)) / ( 4*(1 - lambda)*lambda*m*n) +
    (9*(-2 + 3*n)) / ((1 - lambda)^2*n^2) +
    (6*(26 + 36*m - 54*n - 6*n^2)) / ((1 - lambda)^2*m^2*n^2)

  return(mean.remainderT3)
}


calc.theory.variance.remainderT3.tilde = function(m,n){

  N = n+m
  lambda = m/N

  den.z10_Re = 420*(-1 + lambda)^4 * lambda^4 * m^4 * n^4
  num.z10_Re =
    43740*(4 + 3*m)^2*n^4 +
    4860*lambda*(4 + 3 *m) *n^3 *(21 + 9* m^2 + 7* m*n - 158 *n -
                                    10 *m* (5 + 9 *n)) +
    lambda^4 *(1215 *m^4 *(4 - 3 *n)^2 +
                 108* m^3 *(1568 + 264* n - 1775* n^2 + 1230* n^3) +
                 72* m* (72800 - 14* (12544 + 135* m*n)* n +
                           7* (-11371 + 27 *m*n)* n^2 + (63375 - 1449* m*n)* n^3 + 24924 *n^4) +
                 6* m^2* (625184 + 84* (-4064 + 45* m*n)* n - 15* (38875 + 189* m*n)* n^2 +
                            18810* n^3 + 67464* n^4) +
                 4 *(473200 - 24570* (83 + m*n)* n +
                       21* (50983 + 2916* m*n + 81* m*n^2)* n^2 -
                       1134* (-1993 + 31 *m*n)* n^3 + 487440 *n^4)) +
    9* lambda^2* n^2 *(-424172 + 1215* m^4 + 756* m*n^2 + 749952 *n +
                         659760 *n^2 - 540* m^3 *(61 + 18* n) - 1512 *m*n *(-3 + 32 *n) +
                         54* m^2 *(-7967 + 35 *m*n + 1890 *n + 3890* n^2) -
                         24* m* (39655 - 33534* n - 31608* n^2 + 63 *m*n *(7 + 20* n))) -
    18* lambda^3* n *(-405* m^4* (-4 + 3* n) +
                        18* m^3 *(312 - 480* n + 265* n^2) +
                        4* m^2 *(-7252 + 315* m*n - 118359* n + 10074 *n^2 + 18450* n^3) -
                        2* m* (26992 + 416164 *n - 319539* n^2 - 157746 *n^3 +
                                 63 *m*n *(60 + 36* n + 121* n^2)) +
                        4* (-4095 - 94808* n + 189 *m*n^2 *n + 199773 *n^2 + 75210* n^3 -
                              21* m*n* (65 - 189* n + 291* n^2)))
  z10_Re = num.z10_Re/den.z10_Re

  den.z01_Re = 420*(-1 + lambda)^4 * lambda^4 * m^4 * n^4
  num.z01_Re = 43740 *(4 + 3 *m)^2* n^4 +
    4860 *lambda* (4 + 3 *m)* n^3 *(42 - 64* m + 9* m^2 + 14* m*n - 172* n -
                                      90* m* n) +
    lambda^4* (1215* m^4 *(4 - 3 *n)^2 +
                 18 *m^3* (5824 + 6120 *n - 10335* n^2 + 7380 *n^3) +
                 18 *m* (58240 - 56 *(3904 + 135* m*n)* n -
                           14 *(26411 + 81* m*n) *n^2 + (330801 - 10521* m*n) *n^3 +
                           117210 *n^4) +
                 m^2 *(801472 + 840 *(-2213 + 36 *m*n)* n - 3 *(1486549 + 6615* m*n)* n^2 +
                         272502 *n^3 + 404784 *n^4) +
                 4* (94640 - 24570* (19 + m*n)* n +
                       21 *(-28441 + 2916* m*n + 81* m*n^2)* n^2 -
                       1134 *(-1967 + 61* m*n)* n^3 + 583704* n^4)) +
    81 *lambda^2 *n^2 *(-57484 + 135 *m^4 + 84 *m*n^2 + 89712 *n + 85776* n^2 -
                          90* m^3 *(43 + 12* n) - 168* m*n* (-3 + 62 *n) +
                          3 *m^2 *(-19133 + 105* m*n + 6090* n + 7780* n^2) +
                          2* m *(-64554 + 52881* n + 46610* n^2 - 21 *m*n* (38 + 165 *n))) -
    18 *lambda^3* n *(-405* m^4 *(-4 + 3 *n) + 6* m^3* (796 - 1545*n + 795* n^2) +
                        m^2* (-27748 - 562539 *n + 69906* n^2 + 73800 *n^3 +
                                105* m*n* (16 + 3 *n)) -
                        2* m* (26992 + 521290 *n - 395895* n^2 - 180930* n^3 +
                                 63* m*n* (60 + 66* n + 241* n^2)) +
                        12* (-1365 - 39368* n + 63 *m*n^2* n + 77049* n^2 + 30642 *n^3 -
                               7* m*n* (65 - 189* n + 561* n^2)))

  z01_Re = num.z01_Re/den.z01_Re

  variance.remainder.T3 =  9*(z10_Re*N/m+z01_Re*N/n)

  return(variance.remainder.T3)
}


calc.theory.variance.L.T2.tilde = function(m,n){
  N = n+m
  lambda = m/N

  den = 720*m^2*n^2*lambda^2*(1-lambda)^2
  z10 = (16*n^2 + 2*lambda*n*(14*n+32*m-15)+lambda^2*(75+64*m^2-270*n+256*n^2+56*m*n-60*m))/den
  z01 = (16*n^2 + 2*lambda*n*(14*n+32*m-15)+lambda^2*(15+64*m^2-30*n+16*n^2+56*m*n-60*m))/den

  variance =  4*(z10/lambda + z01/(1-lambda))

  return(variance)
}



exactish.T3.tilde = function(n,m){

  N=n+m
  lambda=m/N

  # mu2.exact
  mu2.Teo1.tilde = compute_theta.Tk.tilde(m=m,n=n,k=2)
  mu2.exact.tilde = mu2.Teo1.tilde + N/(choose(n,2)*choose(m,2))*(N-1)*(N+n)/3

  # mu2.T.tilde
  mu2.L.tilde = 2/(3*n*(1-lambda)) + (2*n-1)/(2*n*m*(1-lambda)) + 2/(3*m*lambda)
  mu2.T.tilde = mu2.Teo1.tilde + mu2.L.tilde

  # mu3.T.tilde
  mu3.Teo1.tilde = compute_theta.Tk.tilde(m=m,n=n,k=3)
  mu3.L.tilde = calc.theory.mean.remainderT3.tilde(m=m,n=n)
  mu3.T.tilde = mu3.Teo1.tilde + mu3.L.tilde

  # mu3.tilde.exact
  mu3.tilde.exact = 2*mu2.exact.tilde + mu2.T.tilde + mu3.T.tilde


  # var2.exact
  var2.Teo1.tilde = calc.theory.variance.Tk.tilde(m=m,n=n,k=2)
  var2.remainder.tilde = calc.theory.variance.remainderT2.tilde(m=m,n=n)
  var2.exactish.tilde = var2.Teo1.tilde + var2.remainder.tilde

  # var2.T.tilde
  var2.L.tilde = calc.theory.variance.L.T2.tilde(m=m,n=n)
  var2.T.tilde = var2.Teo1.tilde + var2.L.tilde

  # var3.T.tilde
  var3.Teo1.tilde = calc.theory.variance.Tk.tilde(m=m,n=n,k=3)
  var3.L.tilde = calc.theory.variance.remainderT3.tilde(m=m,n=n)
  var3.T.tilde = var3.Teo1.tilde + var3.L.tilde

  # var3.tilde.exact
  var3.tilde.exactish = 4*var2.exactish.tilde + var2.T.tilde + var3.T.tilde

  return(c("mu3.tilde.exact" = mu3.tilde.exact,
           "var3.tilde.exactish" = var3.tilde.exactish))
}


exactish.T2.tilde = function(n,m){

  N=n+m
  lambda=m/N

  mu2.tilde.exact = compute_theta.Tk.tilde(m,n,k=2)+ N/(choose(n,2)*choose(m,2))*(N-1)*(N+n)/3


  z01_Re = ((m + n)^2*(16*m^4 + 16*n^4 + 30*m^3*(-1 + 2*n) +
                         30*m*n^2*(-1 + 2*n) + m^2*(15 - 60*n + 92*n^2)))/(45*m^4*n^4)
  z10_Re = ((m + n)^2*(16*m^4 + 16*n^4 + 30*m^3*(-1 + 2*n) +
                         30*m*n^2*(-1 + 2*n) + m^2*(75 + 4*n*(-75 + 83*n))))/(45*m^4*n^4)
  variance.remainder =  4*(z10_Re*N/m+z01_Re*N/n)
  var2.tilde.exactish = calc.theory.variance.Tk.tilde(m,n,k=2) + variance.remainder



  return(c("mu2.tilde.exact" = mu2.tilde.exact,
           "var2.tilde.exactish" = var2.tilde.exactish))
}




