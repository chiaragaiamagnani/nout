library(tidyverse)

n <- 1000
m <- 1000
N <- n+m
lambda <- m/N

gen.data <- function(m,n) {
  Z <- rnorm(m+n)
  return(Z)
}


calc.Tk <- function(Z, m, k) {
  N <- length(Z)
  R <- rank(Z)[(m+1):N]-1
  sum(R^k)
}


calc.Tk.tilde <- function(Z,m,k) {
  N <- length(Z)
  n <- N-m
  Tk <- calc.Tk(Z,m,k)
  #   Tkt <- Tk - factorial(k)^2 * sum(((1:n)-1)^k)
  Tkt <- Tk - sum(((1:n)-1)^k)
  Tkt <- N^(k-1) * Tkt / (choose(n,k) * choose(m,k))
  return(Tkt)
}

# Prende in input la potenza h del binomio e restituisce i coefficienti
# della riga (h+1), che corrisponde ai coefficienti della potenza
# h-esima del binomio
pascalTriangle <- function(h) {
  lapply(0:h, function(i) choose(i, 0:i))[[h+1]]
}


calc.theory.mean.k <- function(m,n,k) {
  N <- n+m
  lambda <- m/N
  r_seq <- 0:(k-1)
  # pascalTriangle(k)[1:k] in questo modo escludo l'ultimo coefficiente
  # che non è incluso nel calcolo della media
  mu <- (factorial(k))^2 * sum( pascalTriangle(k)[1:k] / ( lambda^(r_seq) * (1-lambda)^(k-r_seq-1) )) / (k+1)
  return(mu)
}


aus_e3 = function(k){

  out = 1/((2*k+1))*1/(2^(2*k-2))+
    (k-1)^2 * sum( choose(2*k-2, (0:(2*k-2))) * (-1)^(0:(2*k-2)) / (3:(2*k+1))) +
    (k-1)/(2^(k-2)) * sum(choose(2*k-1, (0:(2*k-1))) * (-1)^(0:(2*k-1)) / (2:(2*k+1)))

  return(out)
}


ausss_E3 = function(k){

  f1 = sum(sapply(1:(k-1),
                  function(h) {
                    sum(sapply(0:(2*h), function(s) choose(2*h,s)*(-1)^s /  (2*k-2*h+s+1)  ) ) /
                      ( (factorial(k-h)*factorial(h))^2 )
                  }
  )
  )

  f2 = 1/((2*k+1)*factorial(k)^2)

  f3 = 2*sum(sapply(1:(k-1),
                    function(h) {
                      sum(sapply(0:(k+h), function(s) choose(k+h,s)*(-1)^s / (k-h+s+1))) /
                        ( factorial(k)*factorial(h)*factorial(k-h) )
                    }
  )
  )

  f4 = 2*sum(sapply(1:(k-1),
                    function(h) {
                      sum(sapply((h+1):(k-1),
                                 function(t){
                                   sum(sapply(0:(h+t), function(s) choose(h+t,s)*(-1)^s / (2*k-h-t+s+1))) /
                                     ( factorial(h)*factorial(k-h)*factorial(t)*factorial(k-t))
                                 }
                      )
                      )
                    }
  )
  )

  out = f1+f2+f3+f4

  return(out)
}




Ehh00 = function(k, lambda){

  out <- (factorial(k))^4/((1-lambda)^(2*k-2))*ausss_E3(k)

  return(out)
}





Ehhr0 = function(k, lambda){

  r_seq <- 1:(k-1)

  coeff_r0 <- (factorial(k))^2 * pascalTriangle(k)[2:k] * factorial(k-r_seq) * factorial(k-r_seq-1) * factorial(r_seq)^2 /
    ( k*lambda^(r_seq)*(1-lambda)^(2*k-r_seq-2) )
  E_r0 <- choose(k-1,k-r_seq-1)*choose(k,1)^2*choose(k-1,r_seq) * ausss_E3(k) +
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
  factorial_r = factorial(k-r_seq) * factorial(k-r_seq-1) * factorial(r_seq)^2
  factorial_s = factorial(k-s_seq) * factorial(k-s_seq-1) * factorial(s_seq)^2
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

  E_rs =  weights1_rs / (k+1)^2 + 2 * weights2_rs / (k+1)^2 + weights3_rs * ausss_E3(k)

  out <- sum(coeff_rs*E_rs)

  return(out)
}



calc.theory.variance.k <- function(m,n,theta,k) {

  N <- n+m
  lambda <- m/N

  if(k==1)
    z_10 <- 1/12

  if(k>1){
    Ehh <- Ehh00(k = k, lambda = lambda) + 2*Ehhr0(k = k, lambda = lambda) + Ehhrs(k = k, lambda = lambda)

    # z_10
    z_10 <- Ehh-theta^2
  }

  # generic k variance
  variance <- z_10*k^2/(N*lambda*(1-lambda))

  return(variance)

}




k <- 2
B <- 1000

n <- 1000
m <- 1000
N <- n+m
lambda <- m/N

calc.theory.mean.k(m,n, k=2)
4*(2-lambda)/(3*lambda*(1-lambda))
calc.theory.variance.k(m,n,theta=calc.theory.mean.k(m,n,k=2),k=2)
64/(45*N*lambda^3*(1-lambda)^3)




# OSSERVAZIONE: più k è grande più la velocità di convergenza è lenta,
# nel senso che n,m devono essere più grandi per arrivare a convergenza per k grande

k <- 1
B <- 1000

n <- 1000
m <- 1000
N <- n+m
lambda <- m/N

stats <- sapply(1:B, function(b) {
  Z <- gen.data(m,n)
  calc.Tk.tilde(Z,m,k)
})

mean(stats)
(mean.theory = calc.theory.mean.k(m,n, k=k))
var(stats)
(variance.theory = calc.theory.variance.k(m,n, theta=mean.theory, k=k))

pp <- tibble(Tkt = stats) %>%
  ggplot(aes(x=Tkt)) +
  geom_histogram(aes(y=..density..)) +
  geom_vline(xintercept=mean.theory, color="red") +
  stat_function(fun = dnorm, args = list(mean = mean.theory, sd = sqrt(variance.theory)),
                color = "blue", lwd = 0.8) +
  theme_bw()
pp


k <- 2
B <- 1000

n <- 1000
m <- 1000
N <- n+m
lambda <- m/N

stats <- sapply(1:B, function(b) {
  Z <- gen.data(m,n)
  calc.Tk.tilde(Z,m,k)
})

mean(stats)
(mean.theory = calc.theory.mean.k(m,n, k=k))
var(stats)
(variance.theory = calc.theory.variance.k(m,n, theta=mean.theory, k=k))

pp <- tibble(Tkt = stats) %>%
  ggplot(aes(x=Tkt)) +
  geom_histogram(aes(y=..density..)) +
  geom_vline(xintercept=mean.theory, color="red") +
  stat_function(fun = dnorm, args = list(mean = mean.theory, sd = sqrt(variance.theory)),
                color = "blue", lwd = 0.8) +
  theme_bw()
pp








k <- 3
B <- 5000

n <- 10000
m <- 10000
N <- n+m
lambda <- m/N

stats <- sapply(1:B, function(b) {
  Z <- gen.data(m,n)
  calc.Tk.tilde(Z,m,k)
})

mean(stats)
(mean.theory = calc.theory.mean.k(m,n, k=k))
var(stats)
(variance.theory = calc.theory.variance.k(m,n, theta=mean.theory, k=k))

pp <- tibble(Tkt = stats) %>%
  ggplot(aes(x=Tkt)) +
  geom_histogram(aes(y=..density..)) +
  geom_vline(xintercept=mean.theory, color="red") +
  stat_function(fun = dnorm, args = list(mean = mean.theory, sd = sqrt(variance.theory)),
                color = "blue", lwd = 0.8) +
  theme_bw()
pp



k <- 5
B <- 5000

n <- 10000
m <- 10000
N <- n+m
lambda <- m/N

stats <- sapply(1:B, function(b) {
  Z <- gen.data(m,n)
  calc.Tk.tilde(Z,m,k)
})

mean(stats)
(mean.theory = calc.theory.mean.k(m,n, k=k))
var(stats)
(variance.theory = calc.theory.variance.k(m,n, theta=mean.theory, k=k))

pp <- tibble(Tkt = stats) %>%
  ggplot(aes(x=Tkt)) +
  geom_histogram() +
  geom_vline(xintercept=mean.theory, color="red") +
  stat_function(fun = dnorm, args = list(mean = mean.theory, sd = sqrt(variance.theory)),
                color = "blue", lwd = 0.8) +
  theme_bw()
pp

















