
# ----------------------------------- NOT EXACT SAMPLE DISTRIBUTION --------------------------------


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


m = 200; n = 200; N = n+m
x = runif(n=m, 0, 2)
y = runif(n=n, 0, 2)
z = c(x,y)


B = 1000*N
WMWstat = vector()
rank3stat = vector()
U = matrix(nrow=B, ncol=n)

sim = function(m, n, z){
  N=m+n
  perm = sample(1:N, m)
  Zperm = c(z[perm], z[-perm])
  # R = rank(Zperm) # ranks in the pooled vector
  U = sapply(1:n, function(j) rank(c(Zperm[1:m], Zperm[m+j]))[m+1]-1)
  return(c("rank3stat" = sum(U^2 + U),
         "WMWstat" = sum(U)))
}

res = foreach(b = 1:B, .combine=cbind) %dopar% {
  values = sim(m=m, n=n, z=z)
  return(values)
}

rank3stat = res[1,]
WMWstat = res[2,]
summary(WMWstat)
summary(rank3stat)

plot(table(WMWstat)/B,
     xlab = "Wilcoxon rank sum statistic",
     ylab = "Density", xaxt="n")
curve(dnorm(x, mean=m*n/2, sd=sqrt((m*n*(m+n+1))/12)), col="blue", add=T)
rug(WMWstat[1])
mean( WMWstat >= WMWstat[1])
axis(side=1, at=WMWstat[1], round(WMWstat[1],3) )




# Approximate standard normal
z11 = 4/45*m^4+16/45*m^3+19/45*m^2+2/15*m
z12 = 4/45*m^4+16/45*m^3+19/45*m^2+2/15*m

theta = n*(m+m*(m-1)/3)
variance = n^2*(m*z11+z12/n)

# LMPI statistic for k=3
plot(table(rank3stat)/length(rank3stat),
     xlab = "LMPI rank statistic for Lehmann's alternative with k=3",
     ylab = "Density",
     ylim = c(0,0.00004),
     xaxt="n",
     yaxt="n")
rug(rank3stat[1])
curve(dnorm(x, mean=theta, sd=sqrt(variance)), col="red", add=T)
axis(side=1, at=rank3stat[1], round(rank3stat[1],3) )
axis(side=2, at=c(0,0.00004), c(0,0.00004))

var(rank3stat)
variance
mean(rank3stat)
theta





# ----------------------------------- EXACT SAMPLE DISTRIBUTION --------------------------------



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
curve(dnorm(x, mean=m*n/2, sd=sqrt((m*n*(m+n+1))/12)), col="blue", add=T)
rug(WMWstat[1])
mean( WMWstat >= WMWstat[1])
axis(side=1, at=WMWstat[1], round(WMWstat[1],3) )



# Approximate standard normal
z11 = 4/45*m^4+16/45*m^3+19/45*m^2+2/15*m
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
curve(dnorm(x, mean=theta, sd=sqrt(variance)), col="red", add=T)
axis(side=1, at=rank3stat[1], round(rank3stat[1],3) )
axis(side=2, at=c(0,0.006), c(0,0.006))

var(rank3stat)
variance
mean(rank3stat)
theta



# Approximate standard normal with sample mean and sample variance
z11 = 4/45*m^4+16/45*m^3+19/45*m^2+2/15*m
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
curve(dnorm(x, mean=mean(rank3stat), sd=sqrt(var(rank3stat))), col="red", add=T)
axis(side=1, at=rank3stat[1], round(rank3stat[1],3) )
axis(side=2, at=c(0,0.006), c(0,0.006))

