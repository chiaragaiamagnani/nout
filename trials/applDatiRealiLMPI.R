library(nout)
library(R.matlab)
set.seed(321)

# Initializing parameters
B=10^3
n = 199
l = 199
m = 20
alpha = m/(l+1)
m1s = seq(from=0, to=m, by=1)

data = readMat("G:\\Il mio Drive\\PHD\\Progetto di ricerca\\Conformal Inference Project\\Simulazioni\\7. Applicazioni dati reali\\Dataset shuttle\\shuttle.mat")
dataset = cbind(data$X, data$y); colnames(dataset)[ncol(dataset)] = "y"
in_ind = which(dataset[,ncol(dataset)]==0)
out_ind = which(dataset[,ncol(dataset)]==1)

res = lapply(m1s,
             function(m1) sim_realdata(B=1, in_index=in_ind, out_index=out_ind,
                                       dataset=dataset,
                                       alpha=alpha,l=l, n=n, m=m, m1=m))

# Storing results
store_res = list("mean.discov" = matrix(nrow=length(m1s), ncol = 5),
                 "mean.powerGlobalNull" = matrix(nrow=length(m1s), ncol = 5))
row.names = rep(NA, times=length(m1s))
for(i in 1:length(m1s)){
  row.names[i] = paste("theta =",m1s[i])
}
rownames(store_res$mean.discov) = row.names
colnames(store_res$mean.discov) = c("BH", "StoBH", "Simes", "StoSimes", "WMW")
rownames(store_res$mean.powerGlobalNull) = row.names
colnames(store_res$mean.powerGlobalNull) = c("BH", "StoBH", "Simes", "StoSimes", "WMW")


for(i in 1:length(res)){
  store_res$mean.discov[i,] = res[[i]]$mean.discov
  store_res$mean.powerGlobalNull[i,] = res[[i]]$mean.powerGlobalNull
}

store_res$mean.discov
store_res$mean.powerGlobalNull

plot(x = m1s, y = store_res$mean.discov[,1], col = 1, ylab = "d",
     xlab = expression(theta), ylim=c(0,m), pch=19,
     main = "Mean of the number of discoveries on B replications")
points(x = m1s, y = store_res$mean.discov[,2], col = 2, pch=19)
points(x = m1s, y = store_res$mean.discov[,5], col = 5, pch=19)
legend("bottomleft", pch = 19, col = c(1,2,5),
       legend =c("BH and Simes CT", "StoreyBH and StoreySimes CT", "WMW CT"))


plot(x = m1s, y = store_res$mean.powerGlobalNull[,1], col = 1, ylab = "power",
     xlab = expression(theta), ylim=c(0,m), pch = 19,
     main = "Mean of the power on B replications")
points(x = m1s, y = store_res$mean.powerGlobalNull[,2], col = 2, pch=19)
points(x = m1s, y = store_res$mean.powerGlobalNull[,5], col = 5, pch=19)
legend("topleft", pch = 19, col = c(1,2,5),
       legend = c("BH and Simes CT", "StoreyBH and StoreySimes CT", "WMW CT"))








