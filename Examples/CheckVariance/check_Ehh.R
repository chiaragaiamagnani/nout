

aus_E3 = function(k){

  f1 = sapply(0:(k-1), function(h) {
    den1 = (factorial(k-h)*factorial(h))^2

    f1.v = sapply(0:(2*k-2*h), function(s) {
      choose(2*k-2*h, s)*(-1)^s/(2*h+s+1)
    })
    sum(f1.v) / den1
  })

  f1.out = sum(f1)


  f2 = lapply(0:(k-1), function(h) {
    if((h+1)<(k-1) || (h+1)==(k-1)){

      sapply((h+1):(k-1), function(t) {
        den2 = factorial(k-h)*factorial(h)*factorial(k-t)*factorial(t)

        f2.v = sapply(0:(2*k-h-t), function(s) {
          choose(2*k-h-t, s)*(-1)^s/(h+s+t+1)
        })
        sum(f2.v) / den2
      })
    }
  })

  f2.vv = sapply(f2, function(x) sum(x))

  f2.out = 2*sum(f2.vv)


  out = f1.out+f2.out

  return(out)
}



hh00 = function(k, m, n){

  N = m+n
  lambda = m/N

  X1 = rnorm(1)
  X = c(X1, rnorm(k-1))
  Y = rnorm(k)
  Xprime = c(X1, rnorm(k-1))
  Yprime = rnorm(k)

  c00 = factorial(k)^4/((1-lambda)^(2*k-2)*k^2)

  mm = sapply(1:k, function(i) X<Y[i])
  p = apply(X = mm, MARGIN = 1, FUN = prod)

  mm.prime = sapply(1:k, function(i) Xprime<Yprime[i])
  pprime = apply(X = mm.prime, MARGIN = 1, FUN = prod)

  v = outer(p, pprime, "*")

  hh00.value = c00*sum(v)

  return(hh00.value)
}



theoric.Ehh00 = function(k, m, n){

  N = m+n
  lambda = m/N

  c00 = factorial(k)^4/((1-lambda)^(2*k-2))

  E3 = ausss_E3(k=k)

  Ehh00 = c00*E3

  return(Ehh00)
}






sample.Ehh00 = mean(replicate(n=10000,h00(k=2,m=100,n=100)))
sample.Ehh00

th.Ehh00 = theoric.Ehh00(k=2,m=100,n=100)
th.Ehh00




sample.Ehh00 = mean(replicate(n=10^5,h00(k=4,m=1000000,n=1000000)))
sample.Ehh00

th.Ehh00 = theoric.Ehh00(k=4,m=100,n=100)
th.Ehh00










kernel.h = function(k, X, Y, lambda){

  j.ind = 1:k
  i.ind = 1:k

  # r=0
  c0 = factorial(k)^2/(k*(1-lambda)^(k-1))
  pi = j.ind
  sigma = combn(i.ind, 1)

  h0 = c0 * sapply(1:length(sigma), function(i) sum(X<Y[sigma[i]]))

  # r>0
  r_seq = 1:(k-1)
  cr = factorial(k-r_seq)*factorial(k-r_seq-1)*factorial(r_seq)^2 / (lambda^r_seq * (1-lambda)^(k-r_seq-1))

  pi_r = sapply(1:length(r_seq), function(i) combn(j.ind, k-r_seq[i]))
  sigma_r = lapply(1:length(r_seq), function(i){
    r = r_seq[i]
    sigma = list()
    for(i0 in i.ind){
      i.ind0 = i.ind[-i0]
      comb_i0 = combn(i.ind0, r)
      sigma[[i0]] = rbind(i0,comb_i0)
    }
    return(sigma)
  }
  )



}




aa = combn(i.ind, k-r_seq[1])


Ehh00 = function(k, lambda){

  out <- (factorial(k))^4/((1-lambda)^(2*k-2))*aus_E3(k)

  return(out)
}





Ehhr0 = function(k, lambda){

  r_seq <- 1:(k-1)

  coeff_r0 <- (factorial(k))^2 * pascalTriangle(k)[2:k] * factorial(k-r_seq) * factorial(k-r_seq-1) * factorial(r_seq) /
    ( k*lambda^(r_seq)*(1-lambda)^(2*k-r_seq-2) )
  E_r0 <- choose(k-1,k-r_seq-1)*choose(k,1)^2*choose(k-1,r_seq) * aus_E3(k) +
    choose(k-1,k-r_seq) * choose(k,1)^2 * choose(k-1,r_seq)/(k+1)^2

  out <- sum(coeff_r0*E_r0)

  return(out)
}




Ehhrs = function(k, lambda){

  r_seq <- 1:(k-1)
  s_seq <- 1:(k-1)

  expo_rs = outer(r_seq, s_seq, "+")
  lambdas = ( lambda^(expo_rs)*(1-lambda)^(2*k-expo_rs-2) )

  c_krs = outer(pascalTriangle(k)[2:k], pascalTriangle(k)[2:k], "*")
  factorial_r = factorial(k) * factorial(k-r_seq) * factorial(k-r_seq-1) * factorial(r_seq)^2 / (k*factorial(k-1))
  factorial_s = factorial(k) * factorial(k-s_seq) * factorial(k-s_seq-1) * factorial(s_seq)^2 / (k*factorial(k-1))
  factorial_rs = outer(factorial_r, factorial_s, "*")
  coeff_rs <- c_krs * factorial_rs / lambdas

  weights1_r = choose(k-1,k-r_seq) * choose(k,1) * choose(k-1,r_seq)
  weights1_s = choose(k-1,k-s_seq) * choose(k,1) * choose(k-1,s_seq)
  weights1_rs = outer(weights1_r, weights1_s, "*")

  weights2_r = choose(k-1,k-r_seq) * choose(k,1) * choose(k-1,r_seq)
  weights2_s = choose(k-1, k-s_seq-1) * choose(k,1) * choose(k-1, s_seq)
  weights2_rs = outer(weights2_r, weights2_s, "*")

  weights3_r = choose(k-1,k-r_seq-1) * choose(k,1) * choose(k-1,r_seq)
  weights3_s = choose(k-1, k-s_seq-1) * choose(k,1) * choose(k-1, s_seq)
  weights3_rs = outer(weights3_r, weights3_s, "*")

  E_rs =  weights1_rs / (k+1)^2 + 2 * weights2_rs / (k+1)^2 + weights3_rs * aus_E3(k)

  out <- sum(coeff_rs*E_rs)

  return(out)
}




