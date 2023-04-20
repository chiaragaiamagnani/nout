# Simu without power / type I error

set.seed(321)

B=10^5
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

# discov = as.data.frame(cbind("d_BH"=d_BH, "d_StoBH"=d_StoBH, "d_Simes"=d_Simes,
#                              "d_StoSimes"=d_StoSimes, "d_WMW"=d_WMW))
discov = as.data.frame(cbind("d_BH"=d_BH>0, "d_StoBH"=d_StoBH>0, "d_Simes"=d_Simes>0,
                             "d_StoSimes"=d_StoSimes>0, "d_WMW"=d_WMW>0))
colnames(discov) = c("BH", "BHSto", "CTSim", "CTSimSto", "CTWMW")
boxplot(discov)
(res = apply(discov, MARGIN = 2, FUN = mean))







thetas = c(0, 1, 0.05, 0.01, 0.07, 0.1)
res = lapply(thetas, function(theta) simuLMPI(B=10^3, n=n, l=l, m=m, d = 3, k = 2, theta, alpha = m/(n+1)))
store_res = matrix(ncol=length(thetas), nrow = 5)
col.names = rep(NA, times=length(thetas))
for(i in 1:length(thetas)){
  col.names[i] = paste("theta =",thetas[i])
}
colnames(store_res) = col.names
rownames(store_res) = c("BH", "StoBH", "Simes", "StoSimes", "WMW")

for(i in 1:length(res)){
  store_res[,i] = res[[i]]$results
}

store_res
















