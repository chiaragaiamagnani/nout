library(mvtnorm)
library(nout)
library(isotree)

set.seed(321)

B=10^2
l = 19
n = 19
m = 2
d = 3
alpha = m/(n+1)

train = as.data.frame(rnorm(n=n))

# H0: F vs H1:G where G = (1-theta)*F + theta*F^2
# Test locally uniform most powerful achieved in correspondance of small values of theta
# F = N(0,1)
# draw_from_mixture = function(theta){
#   success = rbinom(1,1,theta)
#   # if TRUE draw from the alternative distribution
#   if(success==T){
#     x = rmvnorm(n=1, mean=rep(0,d))
#     y = rmvnorm(n=1, mean=rep(0,d))
#     draw = max(x,y)
#   }
#   # if FALSE draw from the null distribution
#   else{
#     draw = rmvnorm(n=1, mean=rep(0,d))
#   }
#   return(draw)
# }


draw_from_mixture = function(theta){
  success = rbinom(1,1,theta)
  # if TRUE draw from the alternative distribution
  if(success==T){
    x = rnorm(n=1)
    y = rnorm(n=1)
    draw = max(x,y)
  }
  # if FALSE draw from the null distribution
  else{
    draw = rnorm(n=1)
  }
  return(draw)
}


iffit = isolation.forest(train, ndim=1, ntrees=10, nthreads=1, scoring_metric = "adj_depth",
                         output_score = TRUE)

crit=critWMW(m=m, n=n, alpha=alpha)
d_WMW = rep(0,B)
d_Simes = rep(0,B)
d_StoSimes = rep(0,B)

for(b in 1:B){
  cal1 = rnorm(n=l)
  cal = as.data.frame(cal1, header = F)
  #aus = draw_from_mixture(theta=0.05)
  te = as.data.frame(unlist(replicate(n=m, draw_from_mixture(theta=0.05))))
  #te = t(prova)

  S_cal = predict.isolation_forest(iffit$model, cal, type = "score")
  S_te = predict.isolation_forest(iffit$model, te, type = "score")

  d_WMW[b]=d_mannwhitney(S_X=S_cal, S_Y=S_te, crit=crit)
  d_Simes[b]=d_Simes(S_X=S_cal, S_Y=S_te, alpha=alpha)
  d_StoSimes[b]=d_StoreySimes(S_X=S_cal, S_Y=S_te, alpha=alpha)
}

















