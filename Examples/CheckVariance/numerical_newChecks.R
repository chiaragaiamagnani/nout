library(tidyverse)

n <- 10000
m <- 10000
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



aus_E3 = function(k){

  f1 = sum(sapply(0:(k-1),
                  function(h) {
                    sum(sapply(0:(2*k-2*h), function(s) choose(2*k-2*h,s)*(-1)^s /  (2*h+s+1)  ) ) /
                      ( (factorial(k-h)*factorial(h))^2 )
                  }))

  f2 = 2*sum(sapply(0:(k-1),
                    function(h) {
                      sum(sapply((h+1):(k-1),
                                 function(t){
                                   sum(sapply(0:(2*k-h-t), function(s) choose(2*k-h-t,s)*(-1)^s / (h+t+s+1))) /
                                     ( factorial(h)*factorial(k-h)*factorial(t)*factorial(k-t))
                                 }))
                      }))

  out = f1+f2

  return(out)
}



Ehh00 = function(k, lambda){

  out <- (factorial(k))^4*aus_E3(k)/((1-lambda)^(2*k-2))

  return(out)
}





Ehhr0 = function(k, lambda){

  r_values = 1:(k-1)
  num_r0 = factorial(k)^2*choose(k,r_values)*factorial(k-r_values)*factorial(k-r_values-1)*(factorial(r_values))^2
  den_r0 = k*lambda^r_values*(1-lambda)^(2*k-r_values-2)
  coeff_r0 = num_r0/den_r0

  f1_r0 = choose(k-1,k-r_values-1)*k^2*choose(k-1,r_values)*aus_E3(k=k)
  f2_r0 = choose(k-1,k-r_values)*k^2*choose(k-1,r_values)/(k+1)^2

  a_r0 = coeff_r0*(f1_r0+f2_r0)

  out = sum(a_r0)

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




# Run experiment
B <- 100
n <- 1000
k <- 2
n.list <- c(10,100,1000)
m.list <- c(0.1,0.2,0.5,1,2,5,10)
L1 = length(n.list)
L2 = length(m.list)
L = L1*L2
n.list.out <- rep(0,L)
lambda.list <- rep(0, L)
empirical.mean.list <- rep(0, L)
theoretical.mean.list <- rep(0, L)
empirical.variance.list <- rep(0, L)
theoretical.variance.list <- rep(0, L)

for(id1 in 1:L1) {
    n <- n.list[id1]
    for(id2 in 1:L2) {
        id = (id1-1)*L2+id2
        m <- round(m.list[id2]*n)
        print(c(id,n,m))
        N <- n+m
        lambda <- m/N
        stats <- sapply(1:B, function(b) {
            Z <- gen.data(m,n)
            calc.Tk.tilde(Z,m)
        })
        mean.theory <- calc.theory.mean.k(m,n)
        variance.theory <- calc.theory.variance.k(m,n)
        n.list.out[id] <- n
        empirical.mean.list[id] <- mean(stats)
        theoretical.mean.list[id] <- mean.theory
        empirical.variance.list[id] <- var(stats)
        theoretical.variance.list[id] <- variance.theory
    }
}
df <- tibble(n=n.list.out,
             Lambda=lambda.list,
             Empirical.mean=empirical.mean.list,
             Theoretical.mean=theoretical.mean.list,
             Empirical.variance=empirical.variance.list,
             Theoretical.variance=theoretical.variance.list)

ppv <- df %>%
  gather(Empirical.variance, Theoretical.variance, key="Variance", value="Value") %>%
  ggplot(aes(x=Lambda, y=Value, color=Variance, shape=Variance)) +
  geom_point() +
  geom_line() +
  scale_y_log10() +
  facet_wrap(.~n, labeller = "label_both") +
  theme_bw()
ppv

# ggsave("plot_mean.pdf", ppv, width=6, height=2.5)
