

library(nout)
library(isotree)
library(foreign)
library(foreach)
library(tidyverse)
library(doSNOW)
library(ggplot2)
library(hommel)
library(mvtnorm)
library(multcomp)

set.seed(123)

# Data from Ludbrook and Dudley (1998)
x = c(5.42, 5.86, 6.16, 6.55, 6.8, 7, 7.11)
y = c(6.51, 7.56, 7.61, 7.84, 11.5)
z = c(x,y)
m = length(x); n = length(y); N = length(z)


m = 10; n = 10; N = n+m
x = runif(n=m, 0, 2)
y = runif(n=n, 0, 2)
z = c(x,y)


perms = combn(1:N,m)
B = ncol(perms) # choose(N,m)
WMWstat = vector()
rank3stat = vector()
U = matrix(nrow=B, ncol=n)
for (b in 1:B){
  Zperm = c(z[perms[,b]], z[-perms[,b]])
  # R = rank(Zperm) # ranks in the pooled vector
  U[b,] = sapply(1:n, function(j) rank(c(Zperm[1:m], Zperm[m+j]))[m+1]-1)
  rank3stat[b] = sum(U[b,]^2 + U[b,])
  WMWstat[b] = sum(U[b,])
}

summary(WMWstat)
summary(rank3stat)

plot(table(WMWstat)/B,
     xlab = "Wilcoxon rank sum statistic",
     ylab = "Frequency", xaxt="n")
lines(dnorm(0:max(WMWstat), mean=m*n/2, sd=sqrt((m*n*(m+n+1))/12)), col="red", add=T)
rug(WMWstat[1])
mean( WMWstat >= WMWstat[1])
axis(side=1, at=WMWstat[1], round(WMWstat[1],3) )




# Approximate standard normal
z11 = 4/45*m^4+16/45*m^3+29/90*m^2+13/30*m-1/5
z12 = 4/45*m^4+16/45*m^3+19/45*m^2+2/15*m

theta = n*(m+m*(m-1)/3)
variance = n^2*(m*z11+z12/n)


# LMPI statistic for k=3
plot(table(rank3stat)/length(rank3stat),
     xlab = "LMPI rank statistic for Lehmann's alternative with k=3",
     ylab = "Frequency",
     ylim = c(0,0.006),
     xaxt="n",
     yaxt="n")
rug(rank3stat[1])
curve(dnorm(x, mean=theta, sd=sqrt(variance)), col="blue", add=T)
lines(dnorm(0:max(rank3stat), mean=theta, sd=sqrt(variance)), col="red")
axis(side=1, at=rank3stat[1], round(rank3stat[1],3) )
axis(side=2, at=c(0,0.006), c(0,0.006))

var(rank3stat)
mean(rank3stat)



plot(table(rank3stat)/length(rank3stat),
     xlab = "LMPI rank statistic for Lehmann's alternative with k=3",
     ylab = "Frequency",
     ylim = c(0,0.006),
     xaxt="n",
     yaxt="n")
rug(rank3stat[1])
curve(dnorm(x, mean=theta, sd=sqrt(22540)), col="red", add=T)
axis(side=1, at=rank3stat[1], round(rank3stat[1],3) )
axis(side=2, at=c(0,0.006), c(0,0.006))




plot(rank3stat,dnorm(rank3stat, mean=theta, sd=sqrt(variance)), col="red")
x11()
