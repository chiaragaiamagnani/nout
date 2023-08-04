# Matrix with lb.d and power for all 5 methods (BH, ABH, Simes, ASimes, WMW) conditional to n1 values

library(nout)
library(R.matlab)
library(isotree)
library(farff)
library(tictoc)
library(tidyverse)
library(doSNOW)
library(ggplot2)
library(hommel)

# -------------------------------- Digits Lehmann's Alternative k=2 -----------------------------------

load("~/nout/trials/RealData/PowerStudy/FinalSimu/Digits/Lehmannk2/resDigits0.1k2")

data = readMat("~/nout/trials/RealData/Datasets/Dataset digits/pendigits.mat")
dataset = cbind(data$X, data$y); colnames(dataset)[ncol(dataset)] = "y"
in_ind = which(dataset[,ncol(dataset)]==0)
out_ind = which(dataset[,ncol(dataset)]==1)

n = 20
n1s = seq(from=0, to=n, by=1)

d_BH = vector()
d_StoBH = vector()
d_Sim = vector()
d_StoSimes = vector()
d_WMW = vector()

pow.rejGlob_BH = vector()
pow.rejGlob_StoBH = vector()
pow.rejGlob_Sim = vector()
pow.rejGlob_StoSimes = vector()
pow.rejGlob_WMW = vector()

for(j in 1:length(n1s)){
  d_BH[j] = resDigits0.1k2$compact.results[[j]]$mean.discoveries[1]
  d_StoBH[j] = resDigits0.1k2$compact.results[[j]]$mean.discoveries[2]
  d_Sim[j] = resDigits0.1k2$compact.results[[j]]$mean.discoveries[3]
  d_StoSimes[j] = resDigits0.1k2$compact.results[[j]]$mean.discoveries[4]
  d_WMW[j] = resDigits0.1k2$compact.results[[j]]$mean.discoveries[5]

  pow.rejGlob_BH[j] = resDigits0.1k2$compact.results[[j]]$mean.powerGlobalNull[1]
  pow.rejGlob_StoBH[j] = resDigits0.1k2$compact.results[[j]]$mean.powerGlobalNull[2]
  pow.rejGlob_Sim[j] = resDigits0.1k2$compact.results[[j]]$mean.powerGlobalNull[3]
  pow.rejGlob_StoSimes[j] = resDigits0.1k2$compact.results[[j]]$mean.powerGlobalNull[4]
  pow.rejGlob_WMW[j] = resDigits0.1k2$compact.results[[j]]$mean.powerGlobalNull[5]

}


lb.d = matrix(nrow = (n+1), ncol = 5)
rownames(lb.d) = as.character(n1s)
colnames(lb.d) = c("FDR-BH", "FDR-Storey", "CT-Simes", "CT-Storey", "CT-WMW")

lb.d[,1] = d_BH
lb.d[,2] = d_StoBH
lb.d[,3] = d_Sim
lb.d[,4] = d_StoSimes
lb.d[,5] = d_WMW
View(lb.d)

pow.rejGlob = matrix(nrow = (n+1), ncol = 5)
rownames(pow.rejGlob) = as.character(seq(from=0, to=n, by=1))
colnames(pow.rejGlob) = c("FDR-BH", "FDR-Storey", "CT-Simes", "CT-Storey", "CT-WMW")
pow.rejGlob[,1] = pow.rejGlob_BH
pow.rejGlob[,2] = pow.rejGlob_StoBH
pow.rejGlob[,3] = pow.rejGlob_Sim
pow.rejGlob[,4] = pow.rejGlob_StoSimes
pow.rejGlob[,5] = pow.rejGlob_WMW
View(pow.rejGlob)

matrixDigits0.1k2 = list("lb.d.matrix" = lb.d, "pow.rejGlob.matrix" = pow.rejGlob)
save(matrixDigits0.1k2, file = "C:/Users/chiar/Documents/matrixDigits0.1k2")

#load("C:/Users/chiar/Documents/matrixDigits0.1k2")
res = matrixDigits0.1k2

theta = seq(0,1, length.out=51)

cat(paste("(",paste(theta,
                    round(sapply(theta, function(p)
                      sum( dbinom(0:n,size=n,prob=p) * res$pow.rejGlob.matrix[,5])
                    ),4), sep=","),")"))

cat(paste("(",paste(theta,
                    round(sapply(theta, function(p)
                      sum( dbinom(0:n,size=n,prob=p) * res$lb.d.matrix[,5])
                    ),4), sep=","),")"))


# -------------------------------- Digits ----------------------------------------------

load("~/nout/trials/RealData/PowerStudy/FinalSimu/Digits/NaturalOutlierDistribution/resDigits0.1")

data = readMat("~/nout/trials/RealData/Datasets/Dataset digits/pendigits.mat")
dataset = cbind(data$X, data$y); colnames(dataset)[ncol(dataset)] = "y"
in_ind = which(dataset[,ncol(dataset)]==0)
out_ind = which(dataset[,ncol(dataset)]==1)

n = 20
n1s = seq(from=0, to=n, by=1)

n.disc_Sim = vector()
n.disc_StoreySim = vector()
n.disc_WMW = vector()

pow.rejGlob_BH = vector()
pow.rejGlob_StoBH = vector()
pow.rejGlob_Sim = vector()
pow.rejGlob_StoSimes = vector()
pow.rejGlob_WMW = vector()

lb.d_BH = vector()
lb.d_StoBH = vector()
lb.d_Sim = vector()
lb.d_StoSimes = vector()
lb.d_WMW = vector()

for(j in 1:length(n1s)){
  lb.d_BH[j] = resDigits0.1$compact.results[[j]]$mean.lb.d[1]
  lb.d_StoBH[j] = resDigits0.1$compact.results[[j]]$mean.lb.d[2]
  lb.d_Sim[j] = resDigits0.1$compact.results[[j]]$mean.lb.d[3]
  lb.d_StoSimes[j] = resDigits0.1$compact.results[[j]]$mean.lb.d[4]
  lb.d_WMW[j] = resDigits0.1$compact.results[[j]]$mean.lb.d[5]

  n.disc_Sim[j] = resDigits0.1$compact.results[[j]]$mean.n.disc[1]
  n.disc_StoreySim[j] = resDigits0.1$compact.results[[j]]$mean.n.disc[3]
  n.disc_WMW[j] = resDigits0.1$compact.results[[j]]$mean.n.disc[4]

  pow.rejGlob_BH[j] = resDigits0.1$compact.results[[j]]$mean.powerGlobalNull[1]
  pow.rejGlob_StoBH[j] = resDigits0.1$compact.results[[j]]$mean.powerGlobalNull[2]
  pow.rejGlob_Sim[j] = resDigits0.1$compact.results[[j]]$mean.powerGlobalNull[3]
  pow.rejGlob_StoSimes[j] = resDigits0.1$compact.results[[j]]$mean.powerGlobalNull[4]
  pow.rejGlob_WMW[j] = resDigits0.1$compact.results[[j]]$mean.powerGlobalNull[5]

}


lb.d = matrix(nrow = (n+1), ncol = 5)
rownames(lb.d) = as.character(n1s)
colnames(lb.d) = c("FDR-BH", "FDR-Storey", "CT-Simes", "CT-Storey", "CT-WMW")

lb.d[,1] = lb.d_BH
lb.d[,2] = lb.d_StoBH
lb.d[,3] = lb.d_Sim
lb.d[,4] = lb.d_StoSimes
lb.d[,5] = lb.d_WMW
View(lb.d)

pow.rejGlob = matrix(nrow = (n+1), ncol = 5)
rownames(pow.rejGlob) = as.character(seq(from=0, to=n, by=1))
colnames(pow.rejGlob) = c("FDR-BH", "FDR-Storey", "CT-Simes", "CT-Storey", "CT-WMW")
pow.rejGlob[,1] = pow.rejGlob_BH
pow.rejGlob[,2] = pow.rejGlob_StoBH
pow.rejGlob[,3] = pow.rejGlob_Sim
pow.rejGlob[,4] = pow.rejGlob_StoSimes
pow.rejGlob[,5] = pow.rejGlob_WMW
View(pow.rejGlob)

n.disc = matrix(nrow = (n+1), ncol = 3)
rownames(n.disc) = as.character(seq(from=0, to=n, by=1))
colnames(n.disc) = c("CT-Simes", "CT-StoreySimes", "CT-WMW")
n.disc[,1] = n.disc_Sim
n.disc[,2] = n.disc_StoreySim
n.disc[,3] = n.disc_WMW
View(n.disc)

matrixDigits0.1 = list("lb.d.matrix" = lb.d, "pow.rejGlob.matrix" = pow.rejGlob, "n.disc" = n.disc)
save(matrixDigits0.1, file = "~/nout/trials/RealData/PowerStudy/FinalSimu/Digits/NaturalOutlierDistribution/matrixDigits0.1")

#load("~/nout/trials/RealData/PowerStudy/FinalSimu/Digits/NaturalOutlierDistribution/matrixDigits0.1")
res = matrixDigits0.1

theta = length(out_ind)/nrow(dataset)

cat(paste("(",paste(round(theta,4),
                    round(sapply(1:5,
                                 function(nc) sum(dbinom(0:n,size=n,prob=theta) * res$pow.rejGlob.matrix[,nc])),
                          4), sep=","),")"))

cat(paste("(",paste(round(theta,4),
                    round(sapply(1:5,
                                 function(nc) sum( dbinom(0:n,size=n,prob=theta) * res$lb.d.matrix[,nc])),
                          4), sep=","),")"))



cat(paste("(",paste(n1s, matrixDigits0.1$lb.d.matrix[,1], sep=","),")"))
cat(paste("(",paste(n1s, matrixDigits0.1$lb.d.matrix[,2], sep=","),")"))
cat(paste("(",paste(n1s, matrixDigits0.1$lb.d.matrix[,3], sep=","),")"))
cat(paste("(",paste(n1s, matrixDigits0.1$lb.d.matrix[,4], sep=","),")"))
cat(paste("(",paste(n1s, matrixDigits0.1$lb.d.matrix[,5], sep=","),")"))


cat(paste("(",paste(n1s, matrixDigits0.1$pow.rejGlob.matrix[,1], sep=","),")"))
cat(paste("(",paste(n1s, matrixDigits0.1$pow.rejGlob.matrix[,2], sep=","),")"))
cat(paste("(",paste(n1s, matrixDigits0.1$pow.rejGlob.matrix[,3], sep=","),")"))
cat(paste("(",paste(n1s, matrixDigits0.1$pow.rejGlob.matrix[,4], sep=","),")"))
cat(paste("(",paste(n1s, matrixDigits0.1$pow.rejGlob.matrix[,5], sep=","),")"))


cat(paste("(",paste(n1s, matrixDigits0.1$n.disc[,1], sep=","),")"))
cat(paste("(",paste(n1s, matrixDigits0.1$n.disc[,2], sep=","),")"))
cat(paste("(",paste(n1s, matrixDigits0.1$n.disc[,3], sep=","),")"))





# -------------------------------- Credit Card -------------------------------------------------

load("~/nout/trials/RealData/PowerStudy/FinalSimu/CreditCard/resCreditCard0.1v2")

dataset = read_csv("~/nout/trials/RealData/Datasets/Dataset creditcard/creditcard.csv")
in_ind = which(dataset[,ncol(dataset)]==0)
out_ind = which(dataset[,ncol(dataset)]==1)

n = 20
n1s = seq(from=0, to=n, by=1)

n.disc_Sim = vector()
n.disc_StoreySim = vector()
n.disc_WMW = vector()

pow.rejGlob_BH = vector()
pow.rejGlob_StoBH = vector()
pow.rejGlob_Sim = vector()
pow.rejGlob_StoSimes = vector()
pow.rejGlob_WMW = vector()

lb.d_BH = vector()
lb.d_StoBH = vector()
lb.d_Sim = vector()
lb.d_StoSimes = vector()
lb.d_WMW = vector()

for(j in 1:length(n1s)){
  lb.d_BH[j] = resCreditCard0.1v2$compact.results[[j]]$mean.lb.d[1]
  lb.d_StoBH[j] = resCreditCard0.1v2$compact.results[[j]]$mean.lb.d[2]
  lb.d_Sim[j] = resCreditCard0.1v2$compact.results[[j]]$mean.lb.d[3]
  lb.d_StoSimes[j] = resCreditCard0.1v2$compact.results[[j]]$mean.lb.d[4]
  lb.d_WMW[j] = resCreditCard0.1v2$compact.results[[j]]$mean.lb.d[5]

  n.disc_Sim[j] = resCreditCard0.1v2$compact.results[[j]]$mean.n.disc[1]
  n.disc_StoreySim[j] = resCreditCard0.1v2$compact.results[[j]]$mean.n.disc[3]
  n.disc_WMW[j] = resCreditCard0.1v2$compact.results[[j]]$mean.n.disc[4]

  pow.rejGlob_BH[j] = resCreditCard0.1v2$compact.results[[j]]$mean.powerGlobalNull[1]
  pow.rejGlob_StoBH[j] = resCreditCard0.1v2$compact.results[[j]]$mean.powerGlobalNull[2]
  pow.rejGlob_Sim[j] = resCreditCard0.1v2$compact.results[[j]]$mean.powerGlobalNull[3]
  pow.rejGlob_StoSimes[j] = resCreditCard0.1v2$compact.results[[j]]$mean.powerGlobalNull[4]
  pow.rejGlob_WMW[j] = resCreditCard0.1v2$compact.results[[j]]$mean.powerGlobalNull[5]

}

lb.d = matrix(nrow = (n+1), ncol = 5)
rownames(lb.d) = as.character(n1s)
colnames(lb.d) = c("FDR-BH", "FDR-Storey", "CT-Simes", "CT-Storey", "CT-WMW")

lb.d[,1] = lb.d_BH
lb.d[,2] = lb.d_StoBH
lb.d[,3] = lb.d_Sim
lb.d[,4] = lb.d_StoSimes
lb.d[,5] = lb.d_WMW
View(lb.d)

pow.rejGlob = matrix(nrow = (n+1), ncol = 5)
rownames(pow.rejGlob) = as.character(seq(from=0, to=n, by=1))
colnames(pow.rejGlob) = c("FDR-BH", "FDR-Storey", "CT-Simes", "CT-Storey", "CT-WMW")
pow.rejGlob[,1] = pow.rejGlob_BH
pow.rejGlob[,2] = pow.rejGlob_StoBH
pow.rejGlob[,3] = pow.rejGlob_Sim
pow.rejGlob[,4] = pow.rejGlob_StoSimes
pow.rejGlob[,5] = pow.rejGlob_WMW
View(pow.rejGlob)

n.disc = matrix(nrow = (n+1), ncol = 3)
rownames(n.disc) = as.character(seq(from=0, to=n, by=1))
colnames(n.disc) = c( "CT-Simes", "CT-StoreySimes", "CT-WMW")
n.disc[,1] = n.disc_Sim
n.disc[,2] = n.disc_StoreySim
n.disc[,3] = n.disc_WMW
View(n.disc)


matrixCredit0.1 = list("lb.d.matrix" = lb.d, "pow.rejGlob.matrix" = pow.rejGlob, "n.disc" = n.disc)
save(matrixCredit0.1, file = "~/nout/trials/RealData/PowerStudy/FinalSimu/CreditCard/matrixCredit0.1")

#load("~/nout/trials/RealData/PowerStudy/FinalSimu/CreditCard/matrixCredit0.1")
res = matrixCredit0.1

theta = length(out_ind)/nrow(dataset)

cat(paste("(",paste(round(theta,4),
                    round(sapply(1:5,
                                 function(nc) sum(dbinom(0:n,size=n,prob=theta) * res$pow.rejGlob.matrix[,nc])),
                          4), sep=","),")"))

cat(paste("(",paste(round(theta,4),
                    round(sapply(1:5,
                                 function(nc) sum( dbinom(0:n,size=n,prob=theta) * res$lb.d.matrix[,nc])),
                          4), sep=","),")"))



cat(paste("(",paste(n1s, matrixCredit0.1$lb.d.matrix[,1], sep=","),")"))
cat(paste("(",paste(n1s, matrixCredit0.1$lb.d.matrix[,2], sep=","),")"))
cat(paste("(",paste(n1s, matrixCredit0.1$lb.d.matrix[,3], sep=","),")"))
cat(paste("(",paste(n1s, matrixCredit0.1$lb.d.matrix[,4], sep=","),")"))
cat(paste("(",paste(n1s, matrixCredit0.1$lb.d.matrix[,5], sep=","),")"))


cat(paste("(",paste(n1s, matrixCredit0.1$pow.rejGlob.matrix[,1], sep=","),")"))
cat(paste("(",paste(n1s, matrixCredit0.1$pow.rejGlob.matrix[,2], sep=","),")"))
cat(paste("(",paste(n1s, matrixCredit0.1$pow.rejGlob.matrix[,3], sep=","),")"))
cat(paste("(",paste(n1s, matrixCredit0.1$pow.rejGlob.matrix[,4], sep=","),")"))
cat(paste("(",paste(n1s, matrixCredit0.1$pow.rejGlob.matrix[,5], sep=","),")"))


cat(paste("(",paste(n1s, matrixCredit0.1$n.disc[,1], sep=","),")"))
cat(paste("(",paste(n1s, matrixCredit0.1$n.disc[,2], sep=","),")"))
cat(paste("(",paste(n1s, matrixCredit0.1$n.disc[,3], sep=","),")"))





# -------------------------------- Shuttle -------------------------------------------------

load("~/nout/trials/RealData/PowerStudy/FinalSimu/Shuttle/resShuttle0.1v2")

data = readMat("~/nout/trials/RealData/Datasets/Dataset shuttle/shuttle.mat")
dataset = cbind(data$X, data$y); colnames(dataset)[ncol(dataset)] = "y"
in_ind = which(dataset[,ncol(dataset)]==0)
out_ind = which(dataset[,ncol(dataset)]==1)

n = 20
n1s = seq(from=0, to=n, by=1)

n.disc_Sim = vector()
n.disc_StoreySim = vector()
n.disc_WMW = vector()

pow.rejGlob_BH = vector()
pow.rejGlob_StoBH = vector()
pow.rejGlob_Sim = vector()
pow.rejGlob_StoSimes = vector()
pow.rejGlob_WMW = vector()

lb.d_BH = vector()
lb.d_StoBH = vector()
lb.d_Sim = vector()
lb.d_StoSimes = vector()
lb.d_WMW = vector()

for(j in 1:length(n1s)){
  lb.d_BH[j] = resShuttle0.1v2$compact.results[[j]]$mean.lb.d[1]
  lb.d_StoBH[j] = resShuttle0.1v2$compact.results[[j]]$mean.lb.d[2]
  lb.d_Sim[j] = resShuttle0.1v2$compact.results[[j]]$mean.lb.d[3]
  lb.d_StoSimes[j] = resShuttle0.1v2$compact.results[[j]]$mean.lb.d[4]
  lb.d_WMW[j] = resShuttle0.1v2$compact.results[[j]]$mean.lb.d[5]

  n.disc_Sim[j] = resShuttle0.1v2$compact.results[[j]]$mean.n.disc[1]
  n.disc_StoreySim[j] = resShuttle0.1v2$compact.results[[j]]$mean.n.disc[3]
  n.disc_WMW[j] = resShuttle0.1v2$compact.results[[j]]$mean.n.disc[4]

  pow.rejGlob_BH[j] = resShuttle0.1v2$compact.results[[j]]$mean.powerGlobalNull[1]
  pow.rejGlob_StoBH[j] = resShuttle0.1v2$compact.results[[j]]$mean.powerGlobalNull[2]
  pow.rejGlob_Sim[j] = resShuttle0.1v2$compact.results[[j]]$mean.powerGlobalNull[3]
  pow.rejGlob_StoSimes[j] = resShuttle0.1v2$compact.results[[j]]$mean.powerGlobalNull[4]
  pow.rejGlob_WMW[j] = resShuttle0.1v2$compact.results[[j]]$mean.powerGlobalNull[5]

}


lb.d = matrix(nrow = (n+1), ncol = 5)
rownames(lb.d) = as.character(n1s)
colnames(lb.d) = c("FDR-BH", "FDR-Storey", "CT-Simes", "CT-Storey", "CT-WMW")

lb.d[,1] = lb.d_BH
lb.d[,2] = lb.d_StoBH
lb.d[,3] = lb.d_Sim
lb.d[,4] = lb.d_StoSimes
lb.d[,5] = lb.d_WMW
View(lb.d)

pow.rejGlob = matrix(nrow = (n+1), ncol = 5)
rownames(pow.rejGlob) = as.character(seq(from=0, to=n, by=1))
colnames(pow.rejGlob) = c("FDR-BH", "FDR-Storey", "CT-Simes", "CT-Storey", "CT-WMW")
pow.rejGlob[,1] = pow.rejGlob_BH
pow.rejGlob[,2] = pow.rejGlob_StoBH
pow.rejGlob[,3] = pow.rejGlob_Sim
pow.rejGlob[,4] = pow.rejGlob_StoSimes
pow.rejGlob[,5] = pow.rejGlob_WMW
View(pow.rejGlob)

n.disc = matrix(nrow = (n+1), ncol = 3)
rownames(n.disc) = as.character(seq(from=0, to=n, by=1))
colnames(n.disc) = c("CT-Simes", "CT-StoreySimes", "CT-WMW")
n.disc[,1] = n.disc_Sim
n.disc[,2] = n.disc_StoreySim
n.disc[,3] = n.disc_WMW
View(n.disc)

matrixShuttle0.1 = list("lb.d.matrix" = lb.d, "pow.rejGlob.matrix" = pow.rejGlob, "n.disc" = n.disc)
save(matrixShuttle0.1, file = "~/nout/trials/RealData/PowerStudy/FinalSimu/Shuttle/matrixShuttle0.1")

#load("~/nout/trials/RealData/PowerStudy/FinalSimu/Shuttle/matrixShuttle0.1")
res = matrixShuttle0.1

theta = length(out_ind)/nrow(dataset)

cat(paste("(",paste(round(theta,4),
                    round(sapply(1:5,
                                 function(nc) sum(dbinom(0:n,size=n,prob=theta) * res$pow.rejGlob.matrix[,nc])),
                          4), sep=","),")"))

cat(paste("(",paste(round(theta,4),
                    round(sapply(1:5,
                                 function(nc) sum( dbinom(0:n,size=n,prob=theta) * res$lb.d.matrix[,nc])),
                          4), sep=","),")"))


cat(paste("(",paste(n1s, matrixShuttle0.1$lb.d.matrix[,1], sep=","),")"))
cat(paste("(",paste(n1s, matrixShuttle0.1$lb.d.matrix[,2], sep=","),")"))
cat(paste("(",paste(n1s, matrixShuttle0.1$lb.d.matrix[,3], sep=","),")"))
cat(paste("(",paste(n1s, matrixShuttle0.1$lb.d.matrix[,4], sep=","),")"))
cat(paste("(",paste(n1s, matrixShuttle0.1$lb.d.matrix[,5], sep=","),")"))


cat(paste("(",paste(n1s, matrixShuttle0.1$pow.rejGlob.matrix[,1], sep=","),")"))
cat(paste("(",paste(n1s, matrixShuttle0.1$pow.rejGlob.matrix[,2], sep=","),")"))
cat(paste("(",paste(n1s, matrixShuttle0.1$pow.rejGlob.matrix[,3], sep=","),")"))
cat(paste("(",paste(n1s, matrixShuttle0.1$pow.rejGlob.matrix[,4], sep=","),")"))
cat(paste("(",paste(n1s, matrixShuttle0.1$pow.rejGlob.matrix[,5], sep=","),")"))


cat(paste("(",paste(n1s, matrixShuttle0.1$n.disc[,1], sep=","),")"))
cat(paste("(",paste(n1s, matrixShuttle0.1$n.disc[,2], sep=","),")"))
cat(paste("(",paste(n1s, matrixShuttle0.1$n.disc[,3], sep=","),")"))





# ----------------------------------- Mammography -------------------------------------

load("~/nout/trials/RealData/PowerStudy/FinalSimu/Mammography/resMammo0.1")

data = readMat("~/nout/trials/RealData/Datasets/Dataset mammography/mammography.mat")
dataset = cbind(data$X, data$y); colnames(dataset)[ncol(dataset)] = "y"
in_ind = which(dataset[,ncol(dataset)]==0)
out_ind = which(dataset[,ncol(dataset)]==1)

n = 20
n1s = seq(from=0, to=n, by=1)

d_BH = vector()
d_StoBH = vector()
d_Sim = vector()
d_StoSimes = vector()
d_WMW = vector()

pow.rejGlob_BH = vector()
pow.rejGlob_StoBH = vector()
pow.rejGlob_Sim = vector()
pow.rejGlob_StoSimes = vector()
pow.rejGlob_WMW = vector()

for(j in 1:length(n1s)){
  d_BH[j] = resMammo0.1$compact.results[[j]]$mean.discoveries[1]
  d_StoBH[j] = resMammo0.1$compact.results[[j]]$mean.discoveries[2]
  d_Sim[j] = resMammo0.1$compact.results[[j]]$mean.discoveries[3]
  d_StoSimes[j] = resMammo0.1$compact.results[[j]]$mean.discoveries[4]
  d_WMW[j] = resMammo0.1$compact.results[[j]]$mean.discoveries[5]

  pow.rejGlob_BH[j] = resMammo0.1$compact.results[[j]]$mean.powerGlobalNull[1]
  pow.rejGlob_StoBH[j] = resMammo0.1$compact.results[[j]]$mean.powerGlobalNull[2]
  pow.rejGlob_Sim[j] = resMammo0.1$compact.results[[j]]$mean.powerGlobalNull[3]
  pow.rejGlob_StoSimes[j] = resMammo0.1$compact.results[[j]]$mean.powerGlobalNull[4]
  pow.rejGlob_WMW[j] = resMammo0.1$compact.results[[j]]$mean.powerGlobalNull[5]

}


lb.d = matrix(nrow = (n+1), ncol = 5)
rownames(lb.d) = as.character(n1s)
colnames(lb.d) = c("FDR-BH", "FDR-Storey", "CT-Simes", "CT-Storey", "CT-WMW")

lb.d[,1] = d_BH
lb.d[,2] = d_StoBH
lb.d[,3] = d_Sim
lb.d[,4] = d_StoSimes
lb.d[,5] = d_WMW
View(lb.d)

pow.rejGlob = matrix(nrow = (n+1), ncol = 5)
rownames(pow.rejGlob) = as.character(seq(from=0, to=n, by=1))
colnames(pow.rejGlob) = c("FDR-BH", "FDR-Storey", "CT-Simes", "CT-Storey", "CT-WMW")
pow.rejGlob[,1] = pow.rejGlob_BH
pow.rejGlob[,2] = pow.rejGlob_StoBH
pow.rejGlob[,3] = pow.rejGlob_Sim
pow.rejGlob[,4] = pow.rejGlob_StoSimes
pow.rejGlob[,5] = pow.rejGlob_WMW
View(pow.rejGlob)

matrixMammo0.1 = list("lb.d.matrix" = lb.d, "pow.rejGlob.matrix" = pow.rejGlob)
save(matrixMammo0.1, file = "~/nout/trials/RealData/PowerStudy/FinalSimu/Mammography/matrixMammo0.1")

#load("~/nout/trials/RealData/PowerStudy/FinalSimu/Mammography/matrixMammo0.1")
res = matrixMammo0.1

theta = length(out_ind)/nrow(dataset)

cat(paste("(",paste(round(theta,4),
                    round(sapply(1:5,
                                 function(nc) sum(dbinom(0:n,size=n,prob=theta) * res$pow.rejGlob.matrix[,nc])),
                          4), sep=","),")"))

cat(paste("(",paste(round(theta,4),
                    round(sapply(1:5,
                                 function(nc)
                                   sum(dbinom(0:n,size=n,prob=theta)*res$lb.d.matrix[,nc])),
                          4), sep=","),")"))

cat(paste("(",paste(n1s, matrixMammo0.1$lb.d.matrix[,1], sep=","),")"))
cat(paste("(",paste(n1s, matrixMammo0.1$lb.d.matrix[,2], sep=","),")"))
cat(paste("(",paste(n1s, matrixMammo0.1$lb.d.matrix[,3], sep=","),")"))
cat(paste("(",paste(n1s, matrixMammo0.1$lb.d.matrix[,4], sep=","),")"))
cat(paste("(",paste(n1s, matrixMammo0.1$lb.d.matrix[,5], sep=","),")"))


cat(paste("(",paste(n1s, matrixMammo0.1$pow.rejGlob.matrix[,1], sep=","),")"))
cat(paste("(",paste(n1s, matrixMammo0.1$pow.rejGlob.matrix[,2], sep=","),")"))
cat(paste("(",paste(n1s, matrixMammo0.1$pow.rejGlob.matrix[,3], sep=","),")"))
cat(paste("(",paste(n1s, matrixMammo0.1$pow.rejGlob.matrix[,4], sep=","),")"))
cat(paste("(",paste(n1s, matrixMammo0.1$pow.rejGlob.matrix[,5], sep=","),")"))



# ---------------------------------- Covertype -----------------------------------------

load("~/nout/trials/RealData/PowerStudy/FinalSimu/Cover/resCover0.1")

data = readMat("~/nout/trials/RealData/Datasets/Dataset cover type/cover.mat")
dataset = cbind(data$X, data$y); colnames(dataset)[ncol(dataset)] = "y"
in_ind = which(dataset[,ncol(dataset)]==0)
out_ind = which(dataset[,ncol(dataset)]==1)

n = 20
n1s = seq(from=0, to=n, by=1)

d_BH = vector()
d_StoBH = vector()
d_Sim = vector()
d_StoSimes = vector()
d_WMW = vector()

pow.rejGlob_BH = vector()
pow.rejGlob_StoBH = vector()
pow.rejGlob_Sim = vector()
pow.rejGlob_StoSimes = vector()
pow.rejGlob_WMW = vector()

for(j in 1:length(n1s)){
  d_BH[j] = resCover0.1$compact.results[[j]]$mean.discoveries[1]
  d_StoBH[j] = resCover0.1$compact.results[[j]]$mean.discoveries[2]
  d_Sim[j] = resCover0.1$compact.results[[j]]$mean.discoveries[3]
  d_StoSimes[j] = resCover0.1$compact.results[[j]]$mean.discoveries[4]
  d_WMW[j] = resCover0.1$compact.results[[j]]$mean.discoveries[5]

  pow.rejGlob_BH[j] = resCover0.1$compact.results[[j]]$mean.powerGlobalNull[1]
  pow.rejGlob_StoBH[j] = resCover0.1$compact.results[[j]]$mean.powerGlobalNull[2]
  pow.rejGlob_Sim[j] = resCover0.1$compact.results[[j]]$mean.powerGlobalNull[3]
  pow.rejGlob_StoSimes[j] = resCover0.1$compact.results[[j]]$mean.powerGlobalNull[4]
  pow.rejGlob_WMW[j] = resCover0.1$compact.results[[j]]$mean.powerGlobalNull[5]

}


lb.d = matrix(nrow = (n+1), ncol = 5)
rownames(lb.d) = as.character(n1s)
colnames(lb.d) = c("FDR-BH", "FDR-Storey", "CT-Simes", "CT-Storey", "CT-WMW")

lb.d[,1] = d_BH
lb.d[,2] = d_StoBH
lb.d[,3] = d_Sim
lb.d[,4] = d_StoSimes
lb.d[,5] = d_WMW
View(lb.d)

pow.rejGlob = matrix(nrow = (n+1), ncol = 5)
rownames(pow.rejGlob) = as.character(seq(from=0, to=n, by=1))
colnames(pow.rejGlob) = c("FDR-BH", "FDR-Storey", "CT-Simes", "CT-Storey", "CT-WMW")
pow.rejGlob[,1] = pow.rejGlob_BH
pow.rejGlob[,2] = pow.rejGlob_StoBH
pow.rejGlob[,3] = pow.rejGlob_Sim
pow.rejGlob[,4] = pow.rejGlob_StoSimes
pow.rejGlob[,5] = pow.rejGlob_WMW
View(pow.rejGlob)

matrixCover0.1 = list("lb.d.matrix" = lb.d, "pow.rejGlob.matrix" = pow.rejGlob)
save(matrixCover0.1, file = "~/nout/trials/RealData/PowerStudy/FinalSimu/Cover/matrixCover0.1")

#load("~/nout/trials/RealData/PowerStudy/FinalSimu/Cover/matrixCover0.1")
res = matrixCover0.1

theta = length(out_ind)/nrow(dataset)

cat(paste("(",paste(round(theta,4),
                    round(sapply(1:5,
                                 function(nc) sum(dbinom(0:n,size=n,prob=theta) * res$pow.rejGlob.matrix[,nc])),
                          4), sep=","),")"))

cat(paste("(",paste(round(theta,4),
                    round(sapply(1:5,
                                 function(nc) sum( dbinom(0:n,size=n,prob=theta) * res$lb.d.matrix[,nc])),
                          4), sep=","),")"))



cat(paste("(",paste(n1s, matrixCover0.1$lb.d.matrix[,1], sep=","),")"))
cat(paste("(",paste(n1s, matrixCover0.1$lb.d.matrix[,2], sep=","),")"))
cat(paste("(",paste(n1s, matrixCover0.1$lb.d.matrix[,3], sep=","),")"))
cat(paste("(",paste(n1s, matrixCover0.1$lb.d.matrix[,4], sep=","),")"))
cat(paste("(",paste(n1s, matrixCover0.1$lb.d.matrix[,5], sep=","),")"))


cat(paste("(",paste(n1s, matrixCover0.1$pow.rejGlob.matrix[,1], sep=","),")"))
cat(paste("(",paste(n1s, matrixCover0.1$pow.rejGlob.matrix[,2], sep=","),")"))
cat(paste("(",paste(n1s, matrixCover0.1$pow.rejGlob.matrix[,3], sep=","),")"))
cat(paste("(",paste(n1s, matrixCover0.1$pow.rejGlob.matrix[,4], sep=","),")"))
cat(paste("(",paste(n1s, matrixCover0.1$pow.rejGlob.matrix[,5], sep=","),")"))


