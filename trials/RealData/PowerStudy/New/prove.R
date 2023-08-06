rm(list=ls())
if(!is.null(dev.list())) dev.off(dev.list()["RStudioGD"])
cat("\014")

library(doSNOW)
library(nout)
library(R.matlab)
library(tictoc)
library(isotree)
library(readr)

dataset = read_csv("G:\\Il mio Drive\\PHD\\Progetto di ricerca\\Conformal Inference Project\\Simulazioni\\7. Applicazioni dati reali\\Dataset creditcard\\creditcard.csv")
out_ind = which(dataset$Class==1)
in_ind = which(dataset$Class==0)

# Initializing parameters
l = 132158
m = 9999
n = 2000
myalpha = n/(m+1)



simParallel = function(B, dataset, in_indexx, out_indexx, n_te, n_outliers, n_cal,
                       iso.ntree, iso.samplesize, alpha, lambda = 0.5, seed=321){

  set.seed(seed)

  tr_ind =  sample(in_indexx, size = l)
  in_ind2 = setdiff(in_indexx, tr_ind)
  tr = dataset[tr_ind,]
  n_cpus = parallel::detectCores()
  iso.fo = isotree::isolation.forest(tr, ndim=ncol(dataset), ntrees=iso.ntree, sample_size = iso.samplesize,
                                     nthreads=n_cpus,
                                     scoring_metric = "depth", output_score = TRUE)
  isofo.model=iso.fo$model

  mycrit = nout::critWMW(m=n,n=m,alpha=alpha)

  # Create the cluster with ncpus-2 cores:
  cluster <- makeCluster(n_cpus-1, type = "SOCK")
  registerDoSNOW(cluster)
  # export dependencies in cluster
  clusterEvalQ(cluster, {c(library(nout), library(hommel), library(isotree), library(dplyr))})
  clusterExport(cluster, list("single_sim", "single_simOut", "dataset", "isofo.model",
                              "alpha", "in_indexx", "out_indexx",
                              "n_te", "n_outliers", "n_cal", "mycrit", "lambda", "B"),
                envir=environment())
  if(n_outliers==0){
    res = snow::parSapply(cl = cluster,
                          X = 1:B,
                          FUN = single_sim,
                          dataset=dataset, model=isofo.model, in_index=in_indexx,
                          out_index=out_indexx, n_outliers=n_outliers, n_te=n_te,
                          n_cal=n_cal, crit=mycrit, alpha=alpha, lambda = lambda)
  }

  else{
    res = snow::parSapply(cl = cluster,
                          X = 1:B,
                          FUN = single_simOut,
                          dataset=dataset, model=isofo.model, in_index=in_indexx,
                          out_index=out_indexx, n_outliers=n_outliers, n_te=n_te,
                          n_cal=n_cal, crit=mycrit, alpha=alpha, lambda = lambda)
  }


  # Stop cluster on master
  stopCluster(cluster)


  return(res)
}





simParallelOut = function(B, dataset, in_index, out_indexx, n_te, n_outliers, n_cal,
                       iso.ntree, iso.samplesize, alpha, lambda = 0.5, seed=321){

  set.seed(seed)

  tr_ind =  sample(in_index, size = l)
  in_ind2 = setdiff(in_index, tr_ind)
  tr = dataset[tr_ind,]
  n_cpus = parallel::detectCores()
  iso.fo = isotree::isolation.forest(tr, ndim=ncol(dataset), ntrees=iso.ntree, sample_size = iso.samplesize,
                                     nthreads=n_cpus,
                                     scoring_metric = "depth", output_score = TRUE)
  isofo.model=iso.fo$model

  mycrit = nout::critWMW(m=n,n=m,alpha=alpha)

  # Create the cluster with ncpus-2 cores:
  cluster <- makeCluster(n_cpus-1, type = "SOCK")
  registerDoSNOW(cluster)
  # export dependencies in cluster
  clusterEvalQ(cluster, {c(library(nout), library(hommel), library(isotree), library(dplyr))})
  clusterExport(cluster, list( "single_simOut", "dataset", "isofo.model"#,
                              #"alpha", "in_index", "out_indexx",
                              #"n_te", "n_outliers", "n_cal", "mycrit", "lambda", "B"
                              ),
                envir=environment()
                )

  res = foreach(b=1:B) %dopar% single_simOut(b, dataset=dataset, model=isofo.model, in_index=as.numeric(in_indexx),
                                          out_index=as.numeric(out_indexx), n_outliers=n_outliers, n_te=n_te,
                                          n_cal=n_cal, crit=mycrit, alpha=alpha, lambda = lambda)
  # res = snow::parSapply(cl = cluster,
  #                         X = 1:B,
  #                         FUN = single_simOut,
  #                         dataset=dataset, model=isofo.model, in_index=in_index,
  #                         out_index=out_index, n_outliers=n_outliers, n_te=n_te,
  #                         n_cal=n_cal, crit=mycrit, alpha=alpha, lambda = lambda)



  # Stop cluster on master
  stopCluster(cluster)


  return(res)
}





single_sim = function(b, dataset, model, in_index, out_index, n_cal, n_te,
                      n_outliers, crit, alpha, lambda = 0.5){

  N=n_te+n_cal
  in_index3 = sample(in_index, size = N)
  cal_ind = in_index3[1:n_cal]
  te_ind = in_index3[(n_cal+1):N]
  cal = dataset[cal_ind,]
  te = dataset[te_ind,]

  S_cal = predict.isolation_forest(model, cal, type = "score")
  S_te = predict.isolation_forest(model, te, type = "score")

  d_WMW = nout::d_mannwhitney(S_Y=S_te, S_X=S_cal, crit = crit)
  d_Sim = nout::d_Simes(S_X=S_cal, S_Y=S_te, alpha=alpha)
  StoSimes = nout::d_StoreySimes(S_X=S_cal, S_Y=S_te, alpha=alpha, lambda=lambda)
  d_StoSimes = StoSimes$d
  pi.not = StoSimes$pi.not
  d_BH = nout::d_benjhoch(S_X=S_cal, S_Y=S_te, alpha=alpha)
  d_StoBH = nout::d_StoreyBH(S_X=S_cal, S_Y=S_te, alpha=alpha, lambda=lambda)
  uniques = length(unique(c(S_cal, S_te)))

  return(as.data.frame(cbind( "d_BH"=d_BH,
                              "d_StoBH"=d_StoBH,
                              "d_Sim"=d_Sim,
                              "d_StoSimes"=d_StoSimes,
                              "d_WMW"=d_WMW,
                              "uniques"=uniques,
                              "n_outliers"=n_outliers,
                              "pi.not"=pi.not,
                              "alpha"=alpha)))
}



single_simOut = function(b, dataset, model, in_index, out_index, n_cal, n_te,
                      n_outliers, crit, alpha, lambda = 0.5){
  n_inliers = n-n_outliers
  N=n_inliers+n_cal
  in_index3 = sample(in_index, size = N)
  cal_ind = in_index3[1:n_cal]
  tein_ind = in_index3[(n_cal+1):N]
  teout_ind = sample(out_index, size = n_outliers)
  cal = dataset[cal_ind,]
  te_ind = c(tein_ind, teout_ind)
  te = dataset[te_ind,]

  S_cal = predict.isolation_forest(model, cal, type = "score")
  S_te = predict.isolation_forest(model, te, type = "score")

  d_WMW = nout::d_mannwhitney(S_Y=S_te, S_X=S_cal, crit = crit)
  d_Sim = nout::d_Simes(S_X=S_cal, S_Y=S_te, alpha=alpha)
  StoSimes = nout::d_StoreySimes(S_X=S_cal, S_Y=S_te, alpha=alpha, lambda=lambda)
  d_StoSimes = StoSimes$d
  pi.not = StoSimes$pi.not
  d_BH = nout::d_benjhoch(S_X=S_cal, S_Y=S_te, alpha=alpha)
  d_StoBH = nout::d_StoreyBH(S_X=S_cal, S_Y=S_te, alpha=alpha, lambda=lambda)
  uniques = length(unique(c(S_cal, S_te)))

  return(as.data.frame(cbind( "d_BH"=d_BH,
                              "d_StoBH"=d_StoBH,
                              "d_Sim"=d_Sim,
                              "d_StoSimes"=d_StoSimes,
                              "d_WMW"=d_WMW,
                              "uniques"=uniques,
                              "n_outliers"=n_outliers,
                              "pi.not"=pi.not,
                              "alpha"=alpha)))
}




compact_results = function(res){
  resT=as.data.frame(t(res))
  colnames(resT) = c("d_BH", "d_StoBH", "d_Simes", "d_StoSimes", "d_WMW", "uniques", "n_outliers", "pi.not", "alpha")

  discoveries = as.data.frame(cbind("d_BH"=unlist(resT$d_BH),
                                    "d_StoBH"=unlist(resT$d_StoBH),
                                    "d_Simes"=unlist(resT$d_Simes),
                                    "d_StoSimes"=unlist(resT$d_StoSimes),
                                    "d_WMW"=unlist(resT$d_WMW)))
  colnames(discoveries) = c("BH", "BHSto", "CTSim", "CTSimSto", "CTWMW")
  mean.discoveries = apply(discoveries, MARGIN = 2, FUN = mean)

  power.GlobalNull = as.data.frame(discoveries>0)
  mean.powerGlobalNull = apply(power.GlobalNull, MARGIN = 2, FUN = mean)

  return(list("discoveries" = discoveries,
              "mean.discoveries" = mean.discoveries,
              "power.GlobalNull" = power.GlobalNull,
              "mean.powerGlobalNull" = mean.powerGlobalNull,
              "pi.not" = unlist(resT$pi.not),
              "uniques"=unlist(resT$uniques),
              "n1"=unlist(resT$n_outliers),
              "alpha"=unlist(resT$alpha)))
}


B=2
re = simParallelOut(B=B, dataset=dataset, in_index=in_ind, out_index=out_ind, iso.ntree=3, iso.samplesize=256,
                 n_te=n, n_outliers=1, n_cal=m, alpha=myalpha, lambda = 0.5, seed=321)



results = compact_results(re)

boxplot(results$discoveries, main="Digits | Distribution of the number of discoveries")
points(x=1:5, y=results$mean.discoveries, pch=19, col="red")
results$mean.discoveries
results$mean.powerGlobalNull



in_index=in_ind
out_index=out_ind
iso.ntree=3
iso.samplesize=256
n_te=n
n_outliers=1
n_cal=m
alpha=myalpha
lambda = 0.5
seed=321
crit=mycrit



somma = function(a,b) a+b
somma(a=3,b="a")
