
library(mvtnorm)
library(nout)
library(isotree)


# scores_from_mixture = function(k, raw_scores, theta){
#
#   if(theta>1 || theta<0){
#     stop("Error: argument theta should in [0,1] interval")
#   }
#
#   ll = length(raw_scores)
#   if(ll<k){
#     stop("Error: length of raw_scores is smaller than k.")
#   }
#
#   quotient = ll%/%k
#   remainder = ll%%k
#
#   if(remainder != 0){
#     cat("Warning: length of raw_scores is not a multiple of k. Last ",
#         remainder, "elements of raw_scores will not be used.")
#   }
#
#   usable.raw_scores = raw_scores[1:(ll-remainder)]
#   success = replicate(quotient, rbinom(1,1,theta))
#
#   scores = rep(0, times = quotient)
#   outlier = rep(0, times = quotient)
#
#   for(i in 0:(quotient-1)){
#     # if TRUE draw from the alternative distribution
#     if(success[i+1]==T){
#       scores[i+1] = max(usable.raw_scores[(i*k+1):(i*k+k)])
#       outlier[i+1]=T
#     }
#     # if FALSE draw from the null distribution
#     if(success[i+1]==F){
#       scores[i+1] = sample(usable.raw_scores[(i*k+1):(i*k+k)], size=1)
#     }
#   }
#
#   return(list("scores"=scores, "outlier"=outlier))
# }



scores_from_mixture = function(k, raw_scores, theta){

  if(theta>1 || theta<0){
    stop("Error: argument theta should in [0,1] interval")
  }

  ll = length(raw_scores)
  if(ll<k){
    stop("Error: length of raw_scores is smaller than k.")
  }

  quotient = ll%/%k # is m
  remainder = ll%%k

  if(remainder != 0){
    cat("Warning: length of raw_scores is not a multiple of k. Last ",
        remainder, "elements of raw_scores will not be used.")
  }

  usable.raw_scores = raw_scores[1:(ll-remainder)]

  m1 = ifelse((theta*m)%%1!=0, round(theta*m), theta*m)

  scores = rep(0, times = quotient)
  outlier = rep(0, times = quotient)

  for(i in 1:m1){
    scores[i+1] = max(usable.raw_scores[(i*k+1):(i*k+k)])
    outlier[i+1]=T
  }

  for(i in (m1+1):m){
    scores[i+1] = sample(usable.raw_scores[(i*k+1):(i*k+k)], size=1)
  }

  return(list("scores"=scores, "outlier"=outlier))
}


simuLMPI = function(B=10^4, n, l, m, d = 3, k = 2, theta, alpha = m/(n+1)){

  train = mvtnorm::rmvnorm(n=n, mean=rep(0,d))
  iso.fo = isotree::isolation.forest(train, ndim=d, ntrees=10, nthreads=1, scoring_metric = "depth",
                                     output_score = TRUE)

  crit=critWMW(m=m, n=n, alpha=alpha)

  d_WMW = rep(0,B)
  d_Simes = rep(0,B)
  d_StoSimes = rep(0,B)
  d_BH = rep(0,B)
  d_StoBH = rep(0,B)

  for(b in 1:B){
    cal = mvtnorm::rmvnorm(n=l, mean=rep(0,d))
    te = mvtnorm::rmvnorm(n=k*m, mean=rep(0,d))

    S_cal = isotree::predict.isolation_forest(iso.fo$model, cal, type = "score")
    rawS_te = isotree::predict.isolation_forest(iso.fo$model, te, type = "score")
    gen.te.score = scores_from_mixture(k=k, raw_scores=rawS_te, theta=theta)
    S_te = gen.te.score$scores

    d_WMW[b] = d_mannwhitney(S_X=S_cal, S_Y=S_te, crit=crit)
    d_Simes[b] = d_Simes(S_X=S_cal, S_Y=S_te, alpha=alpha)
    d_StoSimes[b] = d_StoreySimes(S_X=S_cal, S_Y=S_te, alpha=alpha)
    d_BH[b] = d_benjhoch(S_X=S_cal, S_Y=S_te, alpha=alpha)
    d_StoBH[b] = d_StoreyBH(S_X=S_cal, S_Y=S_te, alpha=alpha)
  }

  discov = as.data.frame(cbind("d_BH"=d_BH, "d_StoBH"=d_StoBH, "d_Simes"=d_Simes,
                               "d_StoSimes"=d_StoSimes, "d_WMW"=d_WMW))
  colnames(discov) = c("BH", "BHSto", "CTSim", "CTSimSto", "CTWMW")
  res = apply(discov, MARGIN = 2, FUN = mean)

  return(list("results"=res, "discoveries"=discov, "theta"=theta, "alpha"=alpha))
}








library(nout)
B=10^3
n = 19
l = 19
m = 2
d = 3
k = 15
alpha = m/(l+1)
m1s = seq(from=0, to=m, by=1)
thetas = m1s/m












