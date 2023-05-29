rm(list=ls())
if(!is.null(dev.list())) dev.off(dev.list()["RStudioGD"])
cat("\014")

# library(Rcpp)
# Rcpp::sourceCpp("~/nout/src/findDiscSum.cpp")

library(doSNOW)
library(foreach)
library(nout)
library(tictoc)
library(isotree)
library(readr)
library(R.matlab)
library(foreign)


compact_results2 = function(res){
  resT=as.data.frame(t(res))

  discoveries = as.data.frame(cbind("dselection_Sim"=unlist(resT$dselection_Sim),
                                    "dselection_WMW"=unlist(resT$dselection_WMW)))
  mean.discoveries = apply(discoveries, MARGIN = 2, FUN = mean)


  power = as.data.frame(discoveries>0)
  mean.power = apply(power, MARGIN = 2, FUN = mean)

  return(list("discoveries" = discoveries,
              "mean.discoveries" = mean.discoveries,
              "power" = power,
              "mean.power" = mean.power,
              "uniques"=unlist(resT$uniques),
              "n1"=unlist(resT$n1),
              "alpha"=unlist(resT$alpha)))
}




data = readMat("~/nout/trials/RealData/Datasets/Dataset digits/pendigits.mat")
dataset = cbind(data$X, data$y); colnames(dataset)[ncol(dataset)] = "y"
in_ind = which(dataset[,ncol(dataset)]==0)
out_ind = which(dataset[,ncol(dataset)]==1)

# Initializing parameters
set.seed(321)

B=10^3

l = 1683
m = 1674
n = 335
m=49
n=10
myalpha = n/(m+1)

tr_ind = sample(in_ind, size = l)
in_ind2 = setdiff(in_ind, tr_ind)
tr = dataset[tr_ind,]
n_cpus = parallel::detectCores()
iso.fo = isotree::isolation.forest(tr, ndim = ncol(dataset), ntrees = 200, sample_size = 256,
                                   nthreads = n_cpus, scoring_metric = "depth",
                                   output_score = TRUE)
isofo.model = iso.fo$model


n1=0

S = 1:n
s = length(S)
n0 = n - n1
N = n0 + m

cl <- makeCluster(parallel::detectCores())
clusterEvalQ(cl, {library(isotree)})
registerDoSNOW(cl)

res = foreach(b = 1:B, .combine=cbind, .packages = c("nout","Rcpp"), .noexport = c("findDiscSum.cpp")) %do% {

  in_index3 = sample(in_ind, size = N)
  cal_ind = in_index3[1:m]
  te_ind = in_index3[(m + 1):N]
  cal = dataset[cal_ind,]
  te = dataset[te_ind,]
  S_X = isotree::predict.isolation_forest(isofo.model, cal, type = "score")
  S_Y = isotree::predict.isolation_forest(isofo.model, te, type = "score")

  dselection_WMW = nout::dselection_mannwhitney(S_Y = S_Y, S_X = S_X,
                                                alpha=myalpha,
                                                S = S)

  dselection_Sim = nout::dselection_Simes(S_X = S_X, S_Y = S_Y, alpha = myalpha,
                                          selection=S)

  uniques = length(unique(c(S_X, S_Y)))
  return(list("dselection_WMW" = dselection_WMW,
              "dselection_Sim" = dselection_Sim,
              "uniques" = uniques,
              "n1" = n1,
              "alpha" = myalpha))
}

stopCluster(cl)


View(res)

res.final = compact_results2(res)





