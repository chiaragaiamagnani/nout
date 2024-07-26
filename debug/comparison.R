
rm(list = ls())


library(tidyverse)
library(latex2exp)
library(gridExtra)
library(foreach)
library(tictoc)
library(fitdistrplus)

library(nout)
source("~/nout/R/utils_g.R")
source("~/nout/mynotes/PowerComparison/utils_simu_g.R")
source("~/nout/R/utils_higher.R")
source("~/nout/R/utils_global_functions.R")
source("~/nout/debug/utils_g_simu_one.R")
source("~/nout/R/utils_LehmannAlt.R")
source("~/nout/debug/utils_g_simu_repeated.R")



## 1) g estimated only once

B=10000; B_MC=500 # number of repetitions
m = 150
n = 150
theta = 0 # proportion of outliers in the test set
alpha = 0.1

rg_null = function(x){ runif(x) }
rg_alt_list = list(function(x){ rg1(x, rg_null=rg_null) },
                   function(x){ rg2(x, rg_null=rg_null) },
                   function(x){ rg9(x, rg_null=rg_null) })
g_alt_list = list(function(x){ g1(x) },
                  function(x){ g2(x) },
                  function(x){ g9(x) })



X = rg_null(x=150)
X1 = X[1:75]
X2 = X[76:150]

# g_KDE = list()
# g_mono = list()


for(i in 1:length(g_alt_list)){
  Y = rg_null(x=75)
  g_mono = estimate_g(X1=X1, X2=X2, Y=Y, constraint="increasing", ker="uniform")
  g_KDE = estimate_g(X1=X1, X2=X2, Y=Y, ker="uniform")

}


# 10 seconds
tic()
res = sim_pow(B=B, B_MC=B_MC, m=m, n=n, theta=theta,
               rg_null=runif, rg=rg1, g=g1,
               g_hat=g_KDE, g_hat_monotone=g_mono,
               alpha=alpha)
toc()




# Generate the data
n.out = as.double(round(theta*n))
n.in = n-n.out
N = as.double(m+n)

X <- rg_null(m)
Y.in <- rg_null(n.in)
if(n.out>0){
  Y.out <- rg(n.out)
  Y <- c(Y.in, Y.out)
} else {
  Y <- Y.in
}
Z <- c(X,Y)

stats_G_N_oracle =  apply(replicate(B_MC, sapply(X=sort(stats::runif(N)), FUN=g)), 1, mean)

T_oracle <- calc.stat.G(Z=Z, m=m, stats_G_vector=stats_G_N_oracle)
T_WMW <- calc.Tk(Z=Z, m=m, k=1)

crit_oracle <-  as.double(stats::qnorm(alpha, mean=meanG(n=n, stats_G_vector=stats_G_N_oracle),
                                       sd = sqrt(varG(n=n, m=m, stats_G_vector=stats_G_N_oracle)), lower.tail = F))
crit_WMW <-  as.double(compute.critical.value.global(m=m, n=n, alpha=alpha, local.test="wmw", k=1, n_perm=0, B=B, seed=123))


T_oracle*(N+1)/2
T_WMW

as.double(stats::qnorm(alpha, mean=(N+1)/2*meanG(n=n, stats_G_vector=stats_G_N_oracle),
                       sd = (N+1)/2*sqrt(varG(n=n, m=m, stats_G_vector=stats_G_N_oracle)), lower.tail = F))
crit_WMW









## 2) g estimated at each iteration

B=100; B_MC=500 # number of repetitions
m=150
n=150
theta = 0
alpha = 0.1


# 480 seconds
tic()
res_rep = sim_pow_rep(B=B, B_MC=B_MC, m=m, n=n, theta=theta,
                      rg_null=runif, rg=rg1, g=g1,
                      alpha=alpha)
toc()



res
res_rep



