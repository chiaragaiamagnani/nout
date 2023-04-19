# Simu with power / type I error

set.seed(321)

B=10^3
n = 19
l = 19
m = 2
d = 3
k = 2
theta = 0.05
alpha = m/(n+1)


train = mvtnorm::rmvnorm(n=n, mean=rep(0,d))
iso.fo = isotree::isolation.forest(train, ndim=d, ntrees=10, nthreads=1, scoring_metric = "depth",
                         output_score = TRUE)

crit=critWMW(m=m, n=n, alpha=alpha)

d_WMW = rep(0,B)
d_Simes = rep(0,B)
d_StoSimes = rep(0,B)
d_BH = rep(0,B)
d_StoBH = rep(0,B)
typeIerr_WMW = rep(0,B)
pow_WMW = rep(0,B)
typeIerr_Simes = rep(0,B)
pow_Simes = rep(0,B)
typeIerr_StoSimes = rep(0,B)
pow_StoSimes = rep(0,B)
typeIerr_BH = rep(0,B)
pow_BH = rep(0,B)
typeIerr_StoBH = rep(0,B)
pow_StoBH = rep(0,B)


for(b in 1:B){
  cal = mvtnorm::rmvnorm(n=l, mean=rep(0,d))
  te = mvtnorm::rmvnorm(n=k*m, mean=rep(0,d))

  S_cal = isotree::predict.isolation_forest(iso.fo$model, cal, type = "score")
  rawS_te = isotree::predict.isolation_forest(iso.fo$model, te, type = "score")
  gen.te.score = scores_from_mixture(k=k, raw_scores=rawS_te, theta=theta)
  S_te = gen.te.score$scores
  out_ind = which(gen.te.score$outlier==1)
  in_ind = which(gen.te.score$outlier==0)

  d_WMW[b] = d_mannwhitney(S_X=S_cal, S_Y=S_te, crit=crit)
  d_Simes[b] = d_Simes(S_X=S_cal, S_Y=S_te, alpha=alpha)
  d_StoSimes[b] = d_StoreySimes(S_X=S_cal, S_Y=S_te, alpha=alpha)
  d_BH[b] = d_benjhoch(S_X=S_cal, S_Y=S_te, alpha=alpha)
  d_StoBH[b] = d_StoreyBH(S_X=S_cal, S_Y=S_te, alpha=alpha)

  if(length(out_ind)==0){
    pow_WMW[b]=0
    pow_Simes[b]=0
    pow_StoSimes[b]=0
    pow_BH[b]=0
    pow_StoBH[b]=0
  }

  if(length(in_ind)==0){
    typeIerr_WMW[b]=0
    typeIerr_Simes[b]=0
    typeIerr_StoSimes[b]=0
    typeIerr_BH[b]=0
    typeIerr_StoBH[b]=0
  }

  if(length(in_ind)!=0 & length(out_ind)!=0){
    typeIerr_WMW[b] = d_mannwhitney(S_X=S_cal, S_Y=S_te[in_ind], crit=crit)>0
    pow_WMW[b] = d_mannwhitney(S_X=S_cal, S_Y=S_te[out_ind], crit=crit)>0

    typeIerr_Simes[b] = d_Simes(S_X=S_cal, S_Y=S_te[in_ind], alpha=alpha)>0
    pow_Simes[b] = d_Simes(S_X=S_cal, S_Y=S_te[out_ind], crit=crit)>0

    typeIerr_StoSimes[b] = d_StoSimes(S_X=S_cal, S_Y=S_te[in_ind], alpha=alpha)>0
    pow_StoSimes[b] = d_StoSimes(S_X=S_cal, S_Y=S_te[out_ind], crit=crit)>0

    typeIerr_BH[b] = d_benjhoch(S_X=S_cal, S_Y=S_te[in_ind], alpha=alpha)>0
    pow_BH[b] = d_benjhoch(S_X=S_cal, S_Y=S_te[out_ind], crit=crit)>0

    typeIerr_StoBH[b] = d_StoreyBH(S_X=S_cal, S_Y=S_te[in_ind], alpha=alpha)>0
    pow_StoBH[b] = d_StoreyBH(S_X=S_cal, S_Y=S_te[out_ind], alpha=alpha)>0
  }
}


discov = as.data.frame(cbind("d_BH"=d_BH, "d_StoBH"=d_StoBH, "d_Simes"=d_Simes, "d_StoSimes"=d_StoSimes, "d_WMW"=d_WMW))
colnames(discov) = c("BH", "BHSto", "CTSim", "CTSimSto", "CTWMW")
boxplot(discov)

pow = as.data.frame(cbind("pow_BH"=pow_BH, "pow_StoBH"=pow_StoBH, "pow_Simes"=pow_Simes,
                          "pow_StoSimes"=pow_StoSimes, "pow_WMW"=pow_WMW))
colnames(pow) = c("BH", "BHSto", "CTSim", "CTSimSto", "CTWMW")
boxplot(pow)

Ierr = as.data.frame(cbind("typeIerr_BH"=typeIerr_BH, "typeIerr_StoBH"=typeIerr_StoBH,
                               "typeIerr_Simes"=typeIerr_Simes, "typeIerr_StoSimes"=typeIerr_StoSimes,
                               "typeIerr_WMW"=typeIerr_WMW))
colnames(Ierr) = c("BH", "BHSto", "CTSim", "CTSimSto", "CTWMW")
boxplot(Ierr)

(power = apply(pow, MARGIN = 2, FUN = mean))
(typeIerr = apply(Ierr, MARGIN = 2, FUN = mean))



