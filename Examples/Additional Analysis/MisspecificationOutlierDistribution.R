
library(tidyverse)

pathDatasets = file.path("G:", "Il mio Drive", "PHD", "Progetto di ricerca",
                         "RankTestsForOutlierDetection",
                         "BiomJ - RankTestsForOutlierDetection", "R code", "Datasets",
                         fsep="\\")

pathResults = file.path("G:", "Il mio Drive", "PHD", "Progetto di ricerca",
                        "RankTestsForOutlierDetection",
                        "BiomJ - RankTestsForOutlierDetection", "R code", fsep="\\")




TrainingIsoForest = function(l, dataset){

  tr_ind = sample(in_ind, size = l)
  tr = dataset[tr_ind,]
  isofo.model = isotree::isolation.forest(tr, ndim=ncol(dataset), ntrees=10, nthreads=1,
                                          scoring_metric = "depth", output_score = TRUE)$model
  in_index2 = setdiff(in_ind, tr_ind)

  return(list("model"=isofo.model, "inlier_remaining" = in_index2))

}



CompareScoreDistributions = function(k, out_ind, inlier_remaining, isofo.model, dataset){

  n1 = length(out_ind)
  n0 = length(inlier_remaining)%/% k
  inliers = dataset[inlier_remaining,]
  outliers = dataset[out_ind[1:n1],]

  # Outlier scores corresponding to Lehmann alternative
  S_inliers = predict.isolation_forest(isofo.model, inliers, type = "score")
  S_outliers.Lehmann = sapply(1:n1, FUN=function(i) max(S_inliers[(1+k*(i-1)):(i*k)]))

  # Outlier scores corresponding to natural outlier distribution
  S_outliers.natural = predict.isolation_forest(isofo.model, outliers, type = "score")

  return(list("S_outliers.natural" = S_outliers.natural,
              "S_outliers.Lehmann" = S_outliers.Lehmann,
              "k" = k))
}



set.seed(321)
l = 1000

data = readMat(paste0(pathDatasets,"\\pendigits.mat"))
dataset = cbind(data$X, data$y); colnames(dataset)[ncol(dataset)] = "y"
in_ind = which(dataset[,ncol(dataset)]==0)
out_ind = which(dataset[,ncol(dataset)]==1)

cluster <- makeCluster(parallel::detectCores())
registerDoSNOW(cluster)
clusterEvalQ(cluster, {list(library(isotree), library(nout))})
clusterExport(cluster, list("l", "in_ind", "out_ind", "dataset", "alpha"))

modeltrain = TrainingIsoForest(l=l, dataset=dataset)
res.scores.k2 = CompareScoreDistributions(k=2,
                          out_ind=out_ind,
                          inlier_remaining=modeltrain$inlier_remaining,
                          isofo.model=modeltrain$model,
                          dataset=dataset)

res.scores.k3 = CompareScoreDistributions(k=3,
                                          out_ind=out_ind,
                                          inlier_remaining=modeltrain$inlier_remaining,
                                          isofo.model=modeltrain$model,
                                          dataset=dataset)

res.scores.k5 = CompareScoreDistributions(k=5,
                                          out_ind=out_ind,
                                          inlier_remaining=modeltrain$inlier_remaining,
                                          isofo.model=modeltrain$model,
                                          dataset=dataset)

res.scores.k8 = CompareScoreDistributions(k=8,
                                          out_ind=out_ind,
                                          inlier_remaining=modeltrain$inlier_remaining,
                                          isofo.model=modeltrain$model,
                                          dataset=dataset)

res.scores.k10 = CompareScoreDistributions(k=10,
                                          out_ind=out_ind,
                                          inlier_remaining=modeltrain$inlier_remaining,
                                          isofo.model=modeltrain$model,
                                          dataset=dataset)

res.scores.k14 = CompareScoreDistributions(k=14,
                                           out_ind=out_ind,
                                           inlier_remaining=modeltrain$inlier_remaining,
                                           isofo.model=modeltrain$model,
                                           dataset=dataset)
stopCluster(cluster)


# Create a new data frame that includes the distribution type for each score
df <- data.frame(
  score = c(res.scores.k2$S_outliers.natural, res.scores.k2$S_outliers.Lehmann),
  distribution = rep(c("natural", "Lehmann"), each = length(res.scores.k2$S_outliers.natural))
)

# Create the plot
pp.scores.k2 <- ggplot(df, aes(x = score, fill = distribution)) +
  geom_histogram(aes(y = ..density..), alpha = 0.5, position = "identity") +
  scale_fill_manual(values = c("natural" = "red", "Lehmann" = "blue")) +
  labs(title = "Natural distribution of outliers vs Lehmann distribution with k=2",
       x = "Scores", y = "Density") +
  theme_bw()

print(pp.scores.k2)


# Create a new data frame that includes the distribution type for each score
df <- data.frame(
  score = c(res.scores.k3$S_outliers.natural, res.scores.k3$S_outliers.Lehmann),
  distribution = rep(c("natural", "Lehmann"), each = length(res.scores.k3$S_outliers.natural))
)

# Create the plot
pp.scores.k3 <- ggplot(df, aes(x = score, fill = distribution)) +
  geom_histogram(aes(y = ..density..), alpha = 0.5, position = "identity") +
  scale_fill_manual(values = c("natural" = "red", "Lehmann" = "blue")) +
  labs(title = "Natural distribution of outliers vs Lehmann distribution with k=3",
       x = "Scores", y = "Density") +
  theme_bw()

print(pp.scores.k3)



# Create a new data frame that includes the distribution type for each score
df <- data.frame(
  score = c(res.scores.k5$S_outliers.natural, res.scores.k5$S_outliers.Lehmann),
  distribution = rep(c("natural", "Lehmann"), each = length(res.scores.k5$S_outliers.natural))
)

# Create the plot
pp.scores.k5 <- ggplot(df, aes(x = score, fill = distribution)) +
  geom_histogram(aes(y = ..density..), alpha = 0.5, position = "identity") +
  scale_fill_manual(values = c("natural" = "red", "Lehmann" = "blue")) +
  labs(title = "Natural distribution of outliers vs Lehmann distribution with k=5",
       x = "Scores", y = "Density") +
  theme_bw()

print(pp.scores.k5)



# Create a new data frame that includes the distribution type for each score
df <- data.frame(
  score = c(res.scores.k8$S_outliers.natural, res.scores.k8$S_outliers.Lehmann),
  distribution = rep(c("natural", "Lehmann"), each = length(res.scores.k8$S_outliers.natural))
)

# Create the plot
pp.scores.k8 <- ggplot(df, aes(x = score, fill = distribution)) +
  geom_histogram(aes(y = ..density..), alpha = 0.5, position = "identity") +
  scale_fill_manual(values = c("natural" = "red", "Lehmann" = "blue")) +
  labs(title = "Natural distribution of outliers vs Lehmann distribution with k=8",
       x = "Scores", y = "Density") +
  theme_bw()

print(pp.scores.k8)



# Create a new data frame that includes the distribution type for each score
df <- data.frame(
  score = c(res.scores.k10$S_outliers.natural, res.scores.k10$S_outliers.Lehmann),
  distribution = rep(c("natural", "Lehmann"), each = length(res.scores.k10$S_outliers.natural))
)

# Create the plot
pp.scores.k10 <- ggplot(df, aes(x = score, fill = distribution)) +
  geom_histogram(aes(y = ..density..), alpha = 0.5, position = "identity") +
  scale_fill_manual(values = c("natural" = "red", "Lehmann" = "blue")) +
  labs(title = "Natural distribution of outliers vs Lehmann distribution with k=10",
       x = "Scores", y = "Density") +
  theme_bw()

print(pp.scores.k10)





# Create a new data frame that includes the distribution type for each score
df <- data.frame(
  score = c(res.scores.k14$S_outliers.natural, res.scores.k14$S_outliers.Lehmann),
  distribution = rep(c("natural", "Lehmann"), each = length(res.scores.k14$S_outliers.natural))
)

# Create the plot
pp.scores.k14 <- ggplot(df, aes(x = score, fill = distribution)) +
  geom_histogram(aes(y = ..density..), alpha = 0.5, position = "identity") +
  scale_fill_manual(values = c("natural" = "red", "Lehmann" = "blue")) +
  labs(title = "Natural distribution of outliers vs Lehmann distribution with k=14",
       x = "Scores", y = "Density") +
  theme_bw()

print(pp.scores.k14)
