
library(nout)
library(R.matlab)
library(isotree)
library(farff)
library(tictoc)
library(tidyverse)
library(doSNOW)
library(ggplot2)


# -------------------------- DIGITS -------------------------------------
data = readMat("~/nout/trials/RealData/Datasets/Dataset digits/pendigits.mat")
dataset = cbind(data$X, data$y); colnames(dataset)[ncol(dataset)] = "y"
in_ind = which(dataset[,ncol(dataset)]==0)
out_ind = which(dataset[,ncol(dataset)]==1)
load("~/nout/trials/RealData/PowerStudy/FinalSimu/Digits/NaturalOutlierDistribution/matrixDigits0.1")


theta = length(out_ind)/nrow(dataset)
n=20
probsn1 = sapply(1:n, function(k) (choose(n,k)*(1-theta)^(n-k)*theta^k)/(1-(1-theta)^n))
power.n1pos = cbind("power.n1pos_BH" = round(sum(matrixDigits0.1$pow.rejGlob.matrix[-1,"FDR-BH"]*probsn1),4),
                    "power.n1pos_StoreyBH" =  round(sum(matrixDigits0.1$pow.rejGlob.matrix[-1,"FDR-Storey"]*probsn1),4),
                    "power.n1pos_WMW" =  round(sum(matrixDigits0.1$pow.rejGlob.matrix[-1,"CT-WMW"]*probsn1),4))

d.n1pos = cbind("d.n1pos_BH" = round(sum(matrixDigits0.1$lb.d.matrix[-1,"FDR-BH"]*probsn1),4),
                "d.n1pos_StoreyBH" =  round(sum(matrixDigits0.1$lb.d.matrix[-1,"FDR-Storey"]*probsn1),4),
                "d.n1pos_Simes" =  round(sum(matrixDigits0.1$lb.d.matrix[-1,"CT-Simes"]*probsn1),4),
                "d.n1pos_CTStorey" =  round(sum(matrixDigits0.1$lb.d.matrix[-1,"CT-Storey"]*probsn1),4),
                "d.n1pos_WMW" =  round(sum(matrixDigits0.1$lb.d.matrix[-1,"CT-WMW"]*probsn1),4))

power.n1is0 = cbind("power.n1is0_BH" = round(matrixDigits0.1$pow.rejGlob.matrix[1,"FDR-BH"],4),
                    "power.n1is0_StoreyBH" =  round(matrixDigits0.1$pow.rejGlob.matrix[1,"FDR-Storey"],4),
                    "power.n1is0_WMW" =  round(matrixDigits0.1$pow.rejGlob.matrix[1,"CT-WMW"],4))

d.n1is0 = cbind("d.n1is0_BH" = round(matrixDigits0.1$lb.d.matrix[1,"FDR-BH"],4),
                "d.n1is0_StoreyBH" =  round(matrixDigits0.1$lb.d.matrix[1,"FDR-Storey"],4),
                "d.n1is0_Simes" =  round(matrixDigits0.1$lb.d.matrix[1,"CT-Simes"],4),
                "d.n1is0_CTStorey" =  round(matrixDigits0.1$lb.d.matrix[1,"CT-Storey"],4),
                "d.n1is0_WMW" =  round(matrixDigits0.1$lb.d.matrix[1,"CT-WMW"],4))




# -------------------------- CREDIT CARD -------------------------------------
dataset = read_csv("~/nout/trials/RealData/Datasets/Dataset creditcard/creditcard.csv")
in_ind = which(dataset[,ncol(dataset)]==0)
out_ind = which(dataset[,ncol(dataset)]==1)
load("~/nout/trials/RealData/PowerStudy/FinalSimu/CreditCard/matrixCredit0.1")


theta = length(out_ind)/nrow(dataset)
n=20
probsn1 = sapply(1:n, function(k) (choose(n,k)*(1-theta)^(n-k)*theta^k)/(1-(1-theta)^n))
power.n1pos = cbind("power.n1pos_BH" = round(sum(matrixCredit0.1$pow.rejGlob.matrix[-1,"FDR-BH"]*probsn1),4),
                    "power.n1pos_StoreyBH" =  round(sum(matrixCredit0.1$pow.rejGlob.matrix[-1,"FDR-Storey"]*probsn1),4),
                    "power.n1pos_WMW" =  round(sum(matrixCredit0.1$pow.rejGlob.matrix[-1,"CT-WMW"]*probsn1),4))

d.n1pos = cbind("d.n1pos_BH" = round(sum(matrixCredit0.1$lb.d.matrix[-1,"FDR-BH"]*probsn1),4),
                "d.n1pos_StoreyBH" =  round(sum(matrixCredit0.1$lb.d.matrix[-1,"FDR-Storey"]*probsn1),4),
                "d.n1pos_Simes" =  round(sum(matrixCredit0.1$lb.d.matrix[-1,"CT-Simes"]*probsn1),4),
                "d.n1pos_CTStorey" =  round(sum(matrixCredit0.1$lb.d.matrix[-1,"CT-Storey"]*probsn1),4),
                "d.n1pos_WMW" =  round(sum(matrixCredit0.1$lb.d.matrix[-1,"CT-WMW"]*probsn1),4))





# -------------------------- SHUTTLE -------------------------------------
data = readMat("~/nout/trials/RealData/Datasets/Dataset shuttle/shuttle.mat")
dataset = cbind(data$X, data$y); colnames(dataset)[ncol(dataset)] = "y"
in_ind = which(dataset[,ncol(dataset)]==0)
out_ind = which(dataset[,ncol(dataset)]==1)
load("~/nout/trials/RealData/PowerStudy/FinalSimu/Shuttle/matrixShuttle0.1")


theta = length(out_ind)/nrow(dataset)
n=20
probsn1 = sapply(1:n, function(k) (choose(n,k)*(1-theta)^(n-k)*theta^k)/(1-(1-theta)^n))
power.n1pos = cbind("power.n1pos_BH" = round(sum(matrixShuttle0.1$pow.rejGlob.matrix[-1,"FDR-BH"]*probsn1),4),
                    "power.n1pos_StoreyBH" =  round(sum(matrixShuttle0.1$pow.rejGlob.matrix[-1,"FDR-Storey"]*probsn1),4),
                    "power.n1pos_WMW" =  round(sum(matrixShuttle0.1$pow.rejGlob.matrix[-1,"CT-WMW"]*probsn1),4))

d.n1pos = cbind("d.n1pos_BH" = round(sum(matrixShuttle0.1$lb.d.matrix[-1,"FDR-BH"]*probsn1),4),
                "d.n1pos_StoreyBH" =  round(sum(matrixShuttle0.1$lb.d.matrix[-1,"FDR-Storey"]*probsn1),4),
                "d.n1pos_Simes" =  round(sum(matrixShuttle0.1$lb.d.matrix[-1,"CT-Simes"]*probsn1),4),
                "d.n1pos_CTStorey" =  round(sum(matrixShuttle0.1$lb.d.matrix[-1,"CT-Storey"]*probsn1),4),
                "d.n1pos_WMW" =  round(sum(matrixShuttle0.1$lb.d.matrix[-1,"CT-WMW"]*probsn1),4))






# -------------------------- COVER -------------------------------------
data = readMat("~/nout/trials/RealData/Datasets/Dataset cover type/cover.mat")
dataset = cbind(data$X, data$y); colnames(dataset)[ncol(dataset)] = "y"
in_ind = which(dataset[,ncol(dataset)]==0)
out_ind = which(dataset[,ncol(dataset)]==1)
load("~/nout/trials/RealData/PowerStudy/FinalSimu/COver/matrixCover0.1")


theta = length(out_ind)/nrow(dataset)
n=20
probsn1 = sapply(1:n, function(k) (choose(n,k)*(1-theta)^(n-k)*theta^k)/(1-(1-theta)^n))
power.n1pos = cbind("power.n1pos_BH" = round(sum(matrixCover0.1$pow.rejGlob.matrix[-1,"FDR-BH"]*probsn1),4),
                    "power.n1pos_StoreyBH" =  round(sum(matrixCover0.1$pow.rejGlob.matrix[-1,"FDR-Storey"]*probsn1),4),
                    "power.n1pos_WMW" =  round(sum(matrixCover0.1$pow.rejGlob.matrix[-1,"CT-WMW"]*probsn1),4))

d.n1pos = cbind("d.n1pos_BH" = round(sum(matrixCover0.1$lb.d.matrix[-1,"FDR-BH"]*probsn1),4),
                "d.n1pos_StoreyBH" =  round(sum(matrixCover0.1$lb.d.matrix[-1,"FDR-Storey"]*probsn1),4),
                "d.n1pos_Simes" =  round(sum(matrixCover0.1$lb.d.matrix[-1,"CT-Simes"]*probsn1),4),
                "d.n1pos_CTStorey" =  round(sum(matrixCover0.1$lb.d.matrix[-1,"CT-Storey"]*probsn1),4),
                "d.n1pos_WMW" =  round(sum(matrixCover0.1$lb.d.matrix[-1,"CT-WMW"]*probsn1),4))




# -------------------------- Mammography -------------------------------------
data = readMat("~/nout/trials/RealData/Datasets/Dataset mammography/mammography.mat")
dataset = cbind(data$X, data$y); colnames(dataset)[ncol(dataset)] = "y"
in_ind = which(dataset[,ncol(dataset)]==0)
out_ind = which(dataset[,ncol(dataset)]==1)
load("~/nout/trials/RealData/PowerStudy/FinalSimu/Mammography/matrixMammo0.1")


theta = length(out_ind)/nrow(dataset)
n=20
probsn1 = sapply(1:n, function(k) (choose(n,k)*(1-theta)^(n-k)*theta^k)/(1-(1-theta)^n))
power.n1pos = cbind("power.n1pos_BH" = round(sum(matrixMammo0.1$pow.rejGlob.matrix[-1,"FDR-BH"]*probsn1),4),
                    "power.n1pos_StoreyBH" =  round(sum(matrixMammo0.1$pow.rejGlob.matrix[-1,"FDR-Storey"]*probsn1),4),
                    "power.n1pos_WMW" =  round(sum(matrixMammo0.1$pow.rejGlob.matrix[-1,"CT-WMW"]*probsn1),4))

d.n1pos = cbind("d.n1pos_BH" = round(sum(matrixMammo0.1$lb.d.matrix[-1,"FDR-BH"]*probsn1),4),
                "d.n1pos_StoreyBH" =  round(sum(matrixMammo0.1$lb.d.matrix[-1,"FDR-Storey"]*probsn1),4),
                "d.n1pos_Simes" =  round(sum(matrixMammo0.1$lb.d.matrix[-1,"CT-Simes"]*probsn1),4),
                "d.n1pos_CTStorey" =  round(sum(matrixMammo0.1$lb.d.matrix[-1,"CT-Storey"]*probsn1),4),
                "d.n1pos_WMW" =  round(sum(matrixMammo0.1$lb.d.matrix[-1,"CT-WMW"]*probsn1),4))






mean(sapply(1:21, function(k) mean(unlist(resMammo0.1$compact.results[[k]]$uniques))))/(199+20)
#mean(sapply(1:21, function(k) mean(unlist(results[[k]]$uniques))))/(199+20)
