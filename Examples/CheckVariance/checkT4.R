#### --------------------------------- CHECK MEAN --------------------------------- ####

library(tidyverse)

n <- 10000
m <- 20000
N <- n+m
lambda <- m/N

gen.data <- function(m,n) {
  Z <- rnorm(m+n)
  return(Z)
}


calc.T4 <- function(Z, m) {
  N <- length(Z)
  R <- rank(Z)[(m+1):N]-1
  sum(R^3)
}

calc.T4.tilde <- function(Z,m) {
  N <- length(Z)
  n <- N-m
  T4 <- calc.T4(Z,m)
  T4t <- T4 - sum(((1:n)-1)^3)
  T4t <- N^2 * T4t / (choose(n,3) * choose(m,3))
  return(T4t)
}

calc.theory.mean.T4 <- function(m,n) {
  N <- n+m
  lambda <- m/N
  mu <- (9*(lambda^2-3*lambda+3))/((1-lambda)^2*lambda^2)
  return(mu)
}


B <- 100
n <- 1000
n.list <- c(10,100,1000)
m.list <- c(0.1,0.2,0.5,1,2,5,10)

L1 = length(n.list)
L2 = length(m.list)
L = L1*L2
n.list.out <- rep(0,L)
lambda.list <- rep(0, L)
empirical.mean.list <- rep(0, L)
theoretical.mean.list <- rep(0, L)
# empirical.variance.list <- rep(0, L)
# theoretical.variance.list <- rep(0, L)

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
      calc.T4.tilde(Z,m)
    })
    mean.theory <- calc.theory.mean.T4(m,n)
    #variance.theory <- calc.theory.variance(m,n)
    n.list.out[id] <- n
    lambda.list[id] <- lambda
    empirical.mean.list[id] <- mean(stats)
    theoretical.mean.list[id] <- mean.theory
    #empirical.variance.list[id] <- var(stats)
    #theoretical.variance.list[id] <- variance.theory
  }
}
df <- tibble(n=n.list.out,
             Lambda=lambda.list,
             Empirical.mean=empirical.mean.list,
             Theoretical.mean=theoretical.mean.list,
             # Empirical.variance=empirical.variance.list,
             # Theoretical.variance=theoretical.variance.list,
)

pp <- df %>%
  gather(Empirical.mean, Theoretical.mean, key="Mean", value="Value") %>%
  ggplot(aes(x=Lambda, y=Value, color=Mean, shape=Mean)) +
  geom_point() +
  geom_line() +
  scale_y_log10() +
  facet_wrap(.~n, labeller = "label_both") +
  theme_bw()
pp




#### --------------------------------- CHECK VARIANCE --------------------------------- ####










