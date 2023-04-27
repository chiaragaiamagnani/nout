
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



scores_from_mixture = function(k, raw_scores, m, m1){

  if(m1==0)  scores = raw_scores

  if(m1==m){
    scores = sapply(0:(m1-1), function(i) max(raw_scores[(i*k+1):(i*k+k)]))
  }

  if(0<m1 & m1<m){
    scores.out = sapply(0:(m1-1), function(i) max(raw_scores[(i*k+1):(i*k+k)]))
    scores.in = raw_scores[(k*m1+1):length(raw_scores)]
    scores = c(scores.out, scores.in)
  }

  return("scores"=scores)
}






simuLMPI = function(B=10^4, n, l, m, d = 3, k = 2, theta, alpha = m/(l+1)){

  train = mvtnorm::rmvnorm(n=n, mean=rep(0,d))
  iso.fo = isotree::isolation.forest(train, ndim=d, ntrees=10, nthreads=1,
                                     scoring_metric = "depth", output_score = TRUE)

  crit=critWMW(m=m, n=n, alpha=alpha)
  m1 = round(theta*m)
  m0 = m-m1

  d_WMW = rep(0,B)
  d_Simes = rep(0,B)
  d_StoSimes = rep(0,B)
  d_BH = rep(0,B)
  d_StoBH = rep(0,B)

  for(b in 1:B){
    cal = mvtnorm::rmvnorm(n=l, mean=rep(0,d))
    te = mvtnorm::rmvnorm(n=k*m1+m0, mean=rep(0,d))

    S_cal = isotree::predict.isolation_forest(iso.fo$model, cal, type = "score")
    rawS_te = isotree::predict.isolation_forest(iso.fo$model, te, type = "score")
    S_te = scores_from_mixture(k=k, raw_scores=rawS_te, m=m, m1=m1)

    d_WMW[b] = d_mannwhitney(S_X=S_cal, S_Y=S_te, crit=crit)
    d_Simes[b] = d_Simes(S_X=S_cal, S_Y=S_te, alpha=alpha)
    d_StoSimes[b] = d_StoreySimes(S_X=S_cal, S_Y=S_te, alpha=alpha)
    d_BH[b] = d_benjhoch(S_X=S_cal, S_Y=S_te, alpha=alpha)
    d_StoBH[b] = d_StoreyBH(S_X=S_cal, S_Y=S_te, alpha=alpha)
  }

  discov = as.data.frame(cbind("d_BH"=d_BH, "d_StoBH"=d_StoBH, "d_Simes"=d_Simes,
                               "d_StoSimes"=d_StoSimes, "d_WMW"=d_WMW))
  colnames(discov) = c("BH", "BHSto", "CTSim", "CTSimSto", "CTWMW")
  mean.discov = apply(discov, MARGIN = 2, FUN = mean)

  powerGlobalNull = as.data.frame(cbind("d_BH"=d_BH>0, "d_StoBH"=d_StoBH>0, "d_Simes"=d_Simes>0,
                                        "d_StoSimes"=d_StoSimes>0, "d_WMW"=d_WMW>0))
  colnames(powerGlobalNull) = c("BH", "BHSto", "CTSim", "CTSimSto", "CTWMW")
  mean.powerGlobalNull = apply(powerGlobalNull, MARGIN = 2, FUN = mean)

  return(list("discoveries"=discov, "mean.discoveries" = mean.discov,
              "powerGlobalNull"=powerGlobalNull, "mean.powerGlobalNull"=mean.powerGlobalNull,
              "theta"=theta, "alpha"=alpha))
}







sim_realdata = function(B, dataset, m1, m, n, l, in_index, out_index=NULL, alpha=m/(l+1), lambda = 0.5){

  m0=m-m1
  if(m1!=0 & is.null(out_index)){
    stop("Error: arg out_index must be initialized.")
  }

  # if(m!=(m1+m0)){
  #   stop("Error: equation m=m1+m0 must be verified.")
  # }

  if(m1!=0){
    tr_ind = sample(in_index, size = n)
    tr = dataset[tr_ind,]
    iso.fo = isolation.forest(tr, ndim=ncol(dataset), ntrees=10, nthreads=1,
                              scoring_metric = "depth", output_score = TRUE)
    in_index2 = setdiff(in_index, tr_ind)

    crit=critWMW(m=m, n=n, alpha=alpha)

    d_WMW = rep(0,B)
    d_Simes = rep(0,B)
    d_StoSimes = rep(0,B)
    d_BH = rep(0,B)
    d_StoBH = rep(0,B)

    for(b in 1:B){
      cal_ind = sample(in_index2, size = l)
      in_index3 = setdiff(in_index2, cal_ind)
      tein_ind = sample(in_index3, size = m0)
      teout_ind = sample(out_index, size = m1)

      cal = dataset[cal_ind,]
      te = dataset[c(tein_ind, teout_ind),]

      S_cal = predict.isolation_forest(iso.fo$model, cal, type = "score")
      S_te = predict.isolation_forest(iso.fo$model, te, type = "score")

      d_WMW[b] = d_mannwhitney(S_X=S_cal, S_Y=S_te, crit=crit)
      d_Simes[b] = d_Simes(S_X=S_cal, S_Y=S_te, alpha=alpha)
      d_StoSimes[b] = d_StoreySimes(S_X=S_cal, S_Y=S_te, alpha=alpha)
      d_BH[b] = d_benjhoch(S_X=S_cal, S_Y=S_te, alpha=alpha)
      d_StoBH[b] = d_StoreyBH(S_X=S_cal, S_Y=S_te, alpha=alpha)
    }
  }

  else{
    tr_ind = sample(in_index, size = n)
    tr = dataset[tr_ind,]
    iso.fo = isolation.forest(tr, ndim=ncol(dataset), ntrees=10, nthreads=1,
                              scoring_metric = "depth", output_score = TRUE)
    in_index2 = setdiff(in_index, tr_ind)

    crit=critWMW(m=m, n=n, alpha=alpha)

    d_WMW = rep(0,B)
    d_Simes = rep(0,B)
    d_StoSimes = rep(0,B)
    d_BH = rep(0,B)
    d_StoBH = rep(0,B)

    for(b in 1:B){
      cal_ind = sample(in_index2, size = l)
      in_index3 = setdiff(in_index2, cal_ind)
      te_ind = sample(in_index3, size = m0)

      cal = dataset[cal_ind,]
      te = dataset[te_ind,]

      S_cal = predict.isolation_forest(iso.fo$model, cal, type = "score")
      S_te = predict.isolation_forest(iso.fo$model, te, type = "score")

      d_WMW[b] = d_mannwhitney(S_X=S_cal, S_Y=S_te, crit=crit)
      d_Simes[b] = d_Simes(S_X=S_cal, S_Y=S_te, alpha=alpha)
      d_StoSimes[b] = d_StoreySimes(S_X=S_cal, S_Y=S_te, alpha=alpha)
      d_BH[b] = d_benjhoch(S_X=S_cal, S_Y=S_te, alpha=alpha)
      d_StoBH[b] = d_StoreyBH(S_X=S_cal, S_Y=S_te, alpha=alpha)
    }
  }

  discov = as.data.frame(cbind("d_BH"=d_BH, "d_StoBH"=d_StoBH, "d_Simes"=d_Simes,
                               "d_StoSimes"=d_StoSimes, "d_WMW"=d_WMW))
  colnames(discov) = c("BH", "BHSto", "CTSim", "CTSimSto", "CTWMW")
  mean.discov = apply(discov, MARGIN = 2, FUN = mean)

  powerGlobalNull = as.data.frame(cbind("d_BH"=d_BH>0, "d_StoBH"=d_StoBH>0, "d_Simes"=d_Simes>0,
                                        "d_StoSimes"=d_StoSimes>0, "d_WMW"=d_WMW>0))
  colnames(powerGlobalNull) = c("BH", "BHSto", "CTSim", "CTSimSto", "CTWMW")
  mean.powerGlobalNull = apply(powerGlobalNull, MARGIN = 2, FUN = mean)

  return(list("discoveries"=discov, "mean.discoveries" = mean.discov,
              "powerGlobalNull"=powerGlobalNull, "mean.powerGlobalNull"=mean.powerGlobalNull,
              "m1"=m1, "alpha"=alpha))

}





132
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












