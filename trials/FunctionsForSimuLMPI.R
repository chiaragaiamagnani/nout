
library(mvtnorm)
library(nout)
library(isotree)



d_benjhoch = function(S_Y, S_X, alpha = 0.1){
  m = length(S_Y)
  n = length(S_X)
  pval = sapply(1:m, function(i) (1+sum(S_X >= S_Y[i]))/(n+1))
  d =  sum(stats::p.adjust(pval,"BH")<=alpha)
  return(d)
}






d_StoreyBH = function(S_Y, S_X, alpha = 0.1, lambda=0.5){
  m = length(S_Y)
  n = length(S_X)
  pval = sort(sapply(1:m, function(i) (1+sum(S_X >= S_Y[i]))/(n+1)), decreasing=FALSE)
  pi0Sto = sapply(1:m, function(i) (1+sum(pval>lambda))/(m*(1-lambda)))
  d =  sum(stats::p.adjust(pval,"BH")<=alpha/pi0Sto)
  return(d)
}




# version 2
#
# d_StoreyBH = function(S_Y, S_X, alpha = 0.1, lambda=0.5){
#   m = length(S_Y)
#   n = length(S_X)
#   pval = sort(sapply(1:m, function(i) (1+sum(S_X >= S_Y[i]))/(n+1)), decreasing=FALSE)
#   pi0Sto = (1+sum(pval>lambda))/(m*(1-lambda))
#   d =  sum(stats::p.adjust(pval,"BH")<=alpha/pi0Sto)
#   return(d)
# }



# version 2 StoreySimes
#
# d_StoreySimes = function(S_Y, S_X, alpha = 0.1, lambda=0.5){
#   m = length(S_Y)
#   n = length(S_X)
#   pval = sort(sapply(1:m, function(i) (1+sum(S_X >= S_Y[i]))/(n+1)), decreasing=FALSE)
#
#   simes.pval = sapply(1:m, function(i)
#     min(pval[i:m]/seq(from=m-i+1, to=1, by=-1)))
#
#   # Building the levels of the Simes test with Storey estimator
#   pi.not = sapply(1:m, function(i)
#     (1+sum(pval[i:m]>lambda))/((m-i+1)*(1-lambda)))
#   coeff = seq(from = m, to = 1, by = -1)
#   thr = alpha/(coeff*pi.not)
#
#   d = sum(cumsum(simes.pval < thr) == 1:m)
#
#   return(d)
# }



scores_from_mixture = function(k, raw_scores, theta){

  if(theta>1 || theta<0){
    stop("Error: argument theta should in [0,1] interval")
  }

  ll = length(raw_scores)
  if(ll<k){
    stop("Error: length of raw_scores is smaller than k.")
  }

  quotient = ll%/%k
  remainder = ll%%k

  if(remainder != 0){
    cat("Warning: length of raw_scores is not a multiple of k. Last ",
        remainder, "elements of raw_scores will not be used.")
  }

  usable.raw_scores = raw_scores[1:(ll-remainder)]
  success = replicate(quotient, rbinom(1,1,theta))

  scores = rep(0, times = quotient)
  outlier = rep(0, times = quotient)

  for(i in 0:(quotient-1)){
    # if TRUE draw from the alternative distribution
    if(success[i+1]==T){
      scores[i+1] = max(usable.raw_scores[(i*k+1):(i*k+k)])
      outlier[i+1]=T
    }
    # if FALSE draw from the null distribution
    if(success[i+1]==F){
      scores[i+1] = sample(usable.raw_scores[(i*k+1):(i*k+k)], size=1)
    }
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










