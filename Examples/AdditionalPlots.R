# Further analysis on Digits dataset

# ---------------------------------------------------------------------------------- #
# 1. Ties in the score vector?
# ---------------------------------------------------------------------------------- #

load("~/nout/Examples/Digits/Lehmannk3/resDigits0.1k3_1999")

m = 1999
n = 200
B=10^4

S_Z = vector()

for(n1 in 1:(n+1)){
  S_Z[n1] = mean(unlist(resDigits0.1k3_1999$raw.res[[n1]]["uniques",]))/(m+n)
}

S_Z




# ---------------------------------------------------------------------------------- #
# 2. Same plot as in Fig.3 of the arxiv with k=2,3,4,5 and theta = 1% and theta = 5%
# ---------------------------------------------------------------------------------- #

library(tidyverse)
library(doSNOW)
library(nout)

gen.data <- function(m,n) {
  Z <- rnorm(m+n)
  return(Z)
}

gen.scores_Lehmann <- function(m, n, n1, k){
  if(n1==0){
    S_Z = gen.data(m,n)
    S_cal = S_Z[1:m]
    S_te = S_Z[(m+1):length(S_Z)]
  }

  if(n1==n){
    augmented.S_Z = gen.data(m,n*k)
    S_cal = augmented.S_Z[1:m]
    augmented.S_te = augmented.S_Z[(m+1):length(augmented.S_Z)]
    S_te = sapply(1:n1, FUN=function(i) max(augmented.S_te[(1+k*(i-1)):(i*k)]))
  }

  if(0<n1&n1<n){
    augmented.S_Z = gen.data(m=m,n=n-n1+n1*k)
    S_cal = augmented.S_Z[1:m]
    augmented.S_te = augmented.S_Z[(m+1):length(augmented.S_Z)]
    inlier.S_te = augmented.S_te[1:(n-n1)]
    outlier.augmented.S_te = augmented.S_te[(n-n1+1):length(augmented.S_te)]
    outlier.S_te = sapply(1:n1, FUN=function(i) max(outlier.augmented.S_te[(1+k*(i-1)):(i*k)]))
    S_te = c(inlier.S_te, outlier.S_te)
  }

  return(list("S_cal" = S_cal,
              "S_te" = S_te,
              "k" = k,
              "n1" = n1))
}



compute_lb.d = function(B, m, n, n1, k, alpha){

  foreach(b = 1:B, .combine=cbind) %dopar% {

    scores = gen.scores_Lehmann(m, n, n1, k)
    S_cal = scores$S_cal
    S_te = scores$S_te
    d_T3 = nout::d_MannWhitneyk3(S_Y = S_te, S_X = S_cal, alpha=alpha)
    d_T3

    d_WMW = nout::d_MannWhitney(S_Y = S_te, S_X = S_cal, alpha=alpha)
    d_T3 = nout::d_MannWhitneyk3(S_Y = S_te, S_X = S_cal, alpha=alpha)
    d_Sim = nout::d_Simes(S_X = S_cal, S_Y = S_te, alpha = alpha)
    StoSimes = nout::d_StoreySimes(S_X = S_cal, S_Y = S_te, alpha = alpha)
    d_StoSimes = StoSimes$d
    pi.not = StoSimes$pi.not
    d_BH = nout::d_benjhoch(S_X = S_cal, S_Y = S_te, alpha = alpha)
    d_StoBH = nout::d_StoreyBH(S_X = S_cal, S_Y = S_te, alpha = alpha)

    return(list("m" = m,
                "n" = n,
                "k" = k,
                "theta" = theta,
                "n1" = n1,
                "alpha" = alpha,
                "S_cal" = S_cal,
                "S_te" = S_te,
                "d_BH" = d_BH,
                "d_StoBH" = d_StoBH,
                "d_Sim" = d_Sim,
                "d_StoSimes" = d_StoSimes,
                "d_WMW" = d_WMW,
                "d_T3" = d_T3,
                "pi.not" = pi.not))
  }
}

CompareT2T3 = function(B, theta, n, m, k, alpha){

  if(theta==0)
    n1 = 0
  if(theta==1)
    n1 = n
  if(0<theta & theta<1)
    n1 = floor(n*theta)

  n1s = seq(from=0, to=n1, by=1)

  res = lapply( 1:length(n1s), function(j) compute_lb.d(B=B, m=m, theta=theta, n=n,
                                                        n1=n1s[j], k=k, alpha=alpha)  )
  return(res)

}


compact_results = function(res, gridd){

  results = list()

  for(j in 1:nrow(gridd)){
    for(theta in 1:length(thetas)){
      lb.d = as.data.frame(
        cbind("d_BH"=unlist(res[[j]][[theta]]["d_BH",]),
              "d_StoBH"=unlist(res[[j]][[theta]]["d_StoBH",]),
              "d_Sim"=unlist(res[[j]][[theta]]["d_Sim",]),
              "d_StoSimes"=unlist(res[[j]][[theta]]["d_StoSimes",]),
              "d_WMW"=unlist(res[[j]][[theta]]["d_WMW",]),
              "d_T3"=unlist(res[[j]][[theta]]["d_T3",])
        )
      )
    mean.lb.d = apply(lb.d, MARGIN = 2, FUN = mean)

    power.GlobalNull = as.data.frame(lb.d>0)
    mean.powerGlobalNull = apply(power.GlobalNull, MARGIN = 2, FUN = mean)

    results[[j]] = list("k" = res[[j]][[theta]]["k",1],
                        "theta" = res[[j]][[theta]]["theta",1],
                        "n" = res[[j]][[theta]]["n",1],
                        "n1" = res[[j]][[theta]]["n1",1],
                        "m" = res[[j]][[theta]]["m",1],
                        "alpha" = res[[j]][[theta]]["alpha",1],
                        "S_cal" = res[[j]][[theta]]["S_cal",],
                        "S_te" = res[[j]][[theta]]["S_te",],
                        "pi.not" = res[[j]][[theta]]["pi.not",],
                        "lb.d" = lb.d,
                        "mean.lb.d" = mean.lb.d,
                        "power.GlobalNull" = power.GlobalNull,
                        "mean.powerGlobalNull" = mean.powerGlobalNull
      )
    }
  }
  return(results)
}







set.seed(321)

B = 10

m = 1999
n = 200
alpha = n/(m+1)

thetas = c(0.01, 0.05, 0.1)
ks = 2:5

# Create a data frame of all combinations of ks and thetas

gridd <- expand.grid(k = ks, theta = thetas)

cluster <- makeCluster(parallel::detectCores()-1)
registerDoSNOW(cluster)
clusterEvalQ(cluster, {list(library(isotree), library(nout))})
clusterExport(cluster, list("n", "m", "thetas", "ks", "alpha", "gen.data", "CompareT2T3", "gen.scores_Lehmann"))

# Apply the CompareT2T3 function to each row of the combinations data frame
res <- lapply(1:nrow(gridd),
              function(i) CompareT2T3(B=B, k=gridd[i, 'k'], theta=gridd[i, 'theta'], n=n, m=m, alpha=alpha))

stopCluster(cluster)


results = compact_results(res, gridd)

d_BH = vector()
d_StoBH = vector()
d_Sim = vector()
d_StoSimes = vector()
d_WMW = vector()
d_T3 = vector()

pow_BH = vector()
pow_StoBH = vector()
pow_Sim = vector()
pow_StoSimes = vector()
pow_WMW = vector()
pow_T3 = vector()

for(j in 1:length(results)){
  d_BH[j] = results[[j]]$mean.lb.d[1]
  d_StoBH[j] = results[[j]]$mean.lb.d[2]
  d_Sim[j] = results[[j]]$mean.lb.d[3]
  d_StoSimes[j] = results[[j]]$mean.lb.d[4]
  d_WMW[j] = results[[j]]$mean.lb.d[5]
  d_T3[j] = results[[j]]$mean.lb.d[6]

  pow_BH[j] = results[[j]]$mean.powerGlobalNull[1]
  pow_StoBH[j] = results[[j]]$mean.powerGlobalNull[2]
  pow_Sim[j] = results[[j]]$mean.powerGlobalNull[3]
  pow_StoSimes[j] = results[[j]]$mean.powerGlobalNull[4]
  pow_WMW[j] = results[[j]]$mean.powerGlobalNull[5]
  pow_T3[j] = results[[j]]$mean.powerGlobalNull[6]
}

# Table unconditional power
thetas = seq(from = 0, to = 1, by = 0.02)
probsn1 = sapply(thetas,
                 function(theta) sapply(1:n,
                                        function(k) choose(n,k)*(1-theta)^(n-k)*theta^(k)))
colnames(probsn1) = as.character(thetas)
rownames(probsn1) = as.character(1:n)
unconditional.power = cbind("uncond.pow_BH" = apply(pow_BH[-1]*probsn1, MARGIN = 2, sum),
                            "uncond.pow_StoreyBH" = apply(pow_StoBH[-1]*probsn1, MARGIN = 2, sum),
                            "uncond.pow_WMW" = apply(pow_WMW[-1]*probsn1, MARGIN = 2, sum),
                            "uncond.pow_T3" = apply(pow_T3[-1]*probsn1, MARGIN = 2, sum))
print(unconditional.power)

d_BH = vector()
d_StoBH = vector()
d_Sim = vector()
d_StoSimes = vector()
d_WMW = vector()
d_T3 = vector()

pow.rejGlob_BH = vector()
pow.rejGlob_StoBH = vector()
pow.rejGlob_Sim = vector()
pow.rejGlob_StoSimes = vector()
pow.rejGlob_WMW = vector()
pow.rejGlob_T3 = vector()


for(j in 1:length(results)){
  d_BH[j] = resDigits0.1k3$compact.results[[j]]$mean.lb.d[1]
  d_StoBH[j] = resDigits0.1k3$compact.results[[j]]$mean.lb.d[2]
  d_Sim[j] = resDigits0.1k3$compact.results[[j]]$mean.lb.d[3]
  d_StoSimes[j] = resDigits0.1k3$compact.results[[j]]$mean.lb.d[4]
  d_WMW[j] = resDigits0.1k3$compact.results[[j]]$mean.lb.d[5]
  d_T3[j] = resDigits0.1k3$compact.results[[j]]$mean.lb.d[6]

  pow.rejGlob_BH[j] = resDigits0.1k3$compact.results[[j]]$mean.powerGlobalNull[1]
  pow.rejGlob_StoBH[j] = resDigits0.1k3$compact.results[[j]]$mean.powerGlobalNull[2]
  pow.rejGlob_Sim[j] = resDigits0.1k3$compact.results[[j]]$mean.powerGlobalNull[3]
  pow.rejGlob_StoSimes[j] = resDigits0.1k3$compact.results[[j]]$mean.powerGlobalNull[4]
  pow.rejGlob_WMW[j] = resDigits0.1k3$compact.results[[j]]$mean.powerGlobalNull[5]
  pow.rejGlob_T3[j] = resDigits0.1k3$compact.results[[j]]$mean.powerGlobalNull[6]

}

lb.d = matrix(nrow = (n+1), ncol = 6)
rownames(lb.d) = as.character(n1s)
colnames(lb.d) = c("FDR-BH", "FDR-Storey", "CT-Simes",
                   "CT-Storey", "CT-WMW", "CT-T3")

lb.d[,1] = d_BH
lb.d[,2] = d_StoBH
lb.d[,3] = d_Sim
lb.d[,4] = d_StoSimes
lb.d[,5] = d_WMW
lb.d[,6] = d_T3

pow.rejGlob = matrix(nrow = (n+1), ncol = 6)
rownames(pow.rejGlob) = as.character(seq(from=0, to=n, by=1))
colnames(pow.rejGlob) = c("FDR-BH", "FDR-Storey", "CT-Simes",
                          "CT-Storey", "CT-WMW", "CT-T3")
pow.rejGlob[,1] = pow.rejGlob_BH
pow.rejGlob[,2] = pow.rejGlob_StoBH
pow.rejGlob[,3] = pow.rejGlob_Sim
pow.rejGlob[,4] = pow.rejGlob_StoSimes
pow.rejGlob[,5] = pow.rejGlob_WMW
pow.rejGlob[,6] = pow.rejGlob_T3

thetas = seq(0,1, length.out=51)

pow_BH = round(sapply(thetas, function(p)
  sum( dbinom(0:n,size=n,prob=p) * res$pow.rejGlob.matrix[,1])),4)
pow_StoBH = round(sapply(thetas, function(p)
  sum( dbinom(0:n,size=n,prob=p) * res$pow.rejGlob.matrix[,2])),4)
pow_Simes = round(sapply(thetas, function(p)
  sum( dbinom(0:n,size=n,prob=p) * res$pow.rejGlob.matrix[,3])),4)
pow_ASimes = round(sapply(thetas, function(p)
  sum( dbinom(0:n,size=n,prob=p) * res$pow.rejGlob.matrix[,4])),4)
pow_WMW = round(sapply(thetas, function(p)
  sum( dbinom(0:n,size=n,prob=p) * res$pow.rejGlob.matrix[,5])),4)
pow_T3 = round(sapply(thetas, function(p)
  sum( dbinom(0:n,size=n,prob=p) * res$pow.rejGlob.matrix[,6])),4)

lb.d.BH = round(sapply(thetas, function(p)
  sum( dbinom(0:n,size=n,prob=p) * res$lb.d.matrix[,1])),4)
lb.d.StoBH = round(sapply(thetas, function(p)
  sum( dbinom(0:n,size=n,prob=p) * res$lb.d.matrix[,2])),4)
lb.d.Simes = round(sapply(thetas, function(p)
  sum( dbinom(0:n,size=n,prob=p) * res$lb.d.matrix[,3])),4)
lb.d.ASimes = round(sapply(thetas, function(p)
  sum( dbinom(0:n,size=n,prob=p) * res$lb.d.matrix[,4])),4)
lb.d.WMW = round(sapply(thetas, function(p)
  sum( dbinom(0:n,size=n,prob=p) * res$lb.d.matrix[,5])),4)
lb.d.T3 = round(sapply(thetas, function(p)
  sum( dbinom(0:n,size=n,prob=p) * res$lb.d.matrix[,6])),4)


# Plot lower bound d
df <- data.frame(
  x = thetas,
  BH = lb.d.BH,
  StoreyBH = lb.d.StoBH,
  Simes_CT = lb.d.Simes,
  StoreySimes_CT = lb.d.ASimes,
  WMW_CT = lb.d.WMW,
  T3_CT = lb.d.T3
)
df_long <- tidyr::pivot_longer(df, cols = -x, names_to = "group", values_to = "y")

ggplot(df_long, aes(x = x, y = y, color = group)) +
  geom_line(size=1) +
  scale_color_manual(values = c("black","gray","blue", "cyan","purple","red")) +
  labs(x = "Number of outliers n1", y = "Average lower bound d") +
  theme_minimal() +
  theme(legend.title = element_blank())


# Plot power
dfpower <- data.frame(
  x = thetas,
  Simes = pow_BH,
  ASimes = pow_StoBH,
  WMW = pow_WMW,
  T3 = pow_T3
)
df_long_power <- tidyr::pivot_longer(dfpower, cols = -x, names_to = "group", values_to = "y")

ggplot(df_long_power, aes(x = x, y = y, color = group)) +
  geom_line(size=1) +
  scale_color_manual(values = c("blue","black","purple","red")) +
  labs(x = "Number of outliers n1", y = "Power to reject the global null") +
  theme_minimal() +
  theme(legend.title = element_blank())

pow_WMW

pow_T3
