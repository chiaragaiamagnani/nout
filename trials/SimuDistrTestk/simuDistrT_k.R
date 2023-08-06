
# SIMULATE T_k STATISTICS PERMUTATION DITRIBUTION

# Library and functions
library(ggplot2)

compute.Tk = function(k, rank.test){
  unsummed = sapply(0:(k-2), function(h) rank.test+h)
  Tk = sum(apply(unsummed, MARGIN = 1, FUN = prod))
  return(Tk)
}


# Data
x = c(5.42, 5.86, 6.16, 6.55, 6.8, 7, 7.11)
y = c(6.51, 7.56, 7.61, 7.84, 11.5)
z = sort(c(x,y), decreasing = F)

m = length(x); n = length(y); N = length(z)

# Compute the value of the T_k (k=7) statistic on rank.y
rank.y = sapply(1:n, function(i) which(z==y[i]))
compute.Tk(k=7, rank.test = rank.y)

# Compute permutation distribution of Uk statistic
B=10^4



# -------------------------------------- k = 7 ----------------------------------------------
Tk.perm = replicate(B, compute.Tk(k=7, rank.test = sample(1:N, n)))

unique.values = unique(Tk.perm)
freqs = sapply(unique.values, function(i) sum(Tk.perm==i)/length(Tk.perm))

# plot 1
df <- data.frame(mids = unique.values, counts = freqs)
ggplot(df, aes(x = mids)) +
  geom_vline(aes(xintercept = mids, alpha = counts), color = "blue") +
  scale_alpha_continuous(range = c(0.1, 1)) +
  ggtitle("LMPI statistic for k = 7") +
  theme_minimal()

# plot 2
df <- data.frame(unique_values = unique.values, freqs = freqs)
ggplot(df, aes(x = unique_values, yend = freqs)) +
  geom_segment(aes(xend = unique_values, y = 0), color = "blue") +
  scale_y_continuous(limits = c(0, max(df$freqs))) +
  ggtitle("LMPI statistic for k = 7") +
  theme_minimal()


# -------------------------------- k = 2 --------------------------------------------
Tk.perm <- replicate(B, compute.Tk(k = 2, rank.test = sample(1:N, n))) - n * (n + 1) / 2

unique.values = sort(unique(Tk.perm), decreasing = F)

freqs = sapply(unique.values, function(i) sum(Tk.perm==i)/length(Tk.perm))

# plot 1
df <- data.frame(mids = unique.values, counts = freqs)
ggplot(df, aes(x = mids)) +
  geom_vline(aes(xintercept = mids, alpha = counts), color = "blue") +
  scale_alpha_continuous(range = c(0.1, 1)) +
  ggtitle("LMPI statistic for k = 2") +
  theme_minimal()

# plot 2
df <- data.frame(unique_values = unique.values, freqs = freqs)
ggplot(df, aes(x = unique_values, yend = freqs)) +
  geom_segment(aes(xend = unique_values, y = 0), color = "blue") +
  scale_y_continuous(limits = c(0, 0.07)) + #limits = c(0, max(df$freqs))
  ggtitle("LMPI statistic for k = 2") +
  stat_function(fun = dnorm, args = list(mean = m*n/2, sd = sqrt((m*n*(m+n+1))/12)), color = "darkgreen") +
  geom_line(data = data.frame(x = unique.values, y = dwilcox(unique.values, m, n)), aes(x = x, y = y), color = "red", size = 1) +
  theme_minimal()



# ----------------------------------- k = 3 -----------------------------------------------
Tk.perm = replicate(B, compute.Tk(k=3, rank.test = sample(1:N, n))) - 3*n*(n-1) -1

unique.values = unique(Tk.perm)
freqs = sapply(unique.values, function(i) sum(Tk.perm==i)/length(Tk.perm))

# plot 1
df <- data.frame(mids = unique.values, counts = freqs)
ggplot(df, aes(x = mids)) +
  geom_vline(aes(xintercept = mids, alpha = counts), color = "blue") +
  scale_alpha_continuous(range = c(0.1, 1)) +
  ggtitle("LMPI statistic for k = 3") +
  theme_minimal()

# Approximate standard normal
z11 = 4/45*m^4+16/45*m^3+29/90*m^2+13/30*m-1/5
z12 = 4/45*m^4+16/45*m^3+19/45*m^2+2/15*m

theta = (m+m*(m-1)/3)
variance = (m*z11+z12/n)/n^2

# Approximate standard normal
z11 = 4/45*N^4+16/45*N^3+29/90*N^2+13/30*N-1/5
z12 = 4/45*N^4+16/45*N^3+19/45*N^2+2/15*N

theta = (N+N*(N-1)/3)
variance = (N*z11+z12/n)/n^2
#quantiles = stats::qnorm(alpha, mean=theta, sd = sqrt(variance), lower.tail = F)

# plot 2
df <- data.frame(unique_values = unique.values, freqs = freqs)
ggplot(df, aes(x = unique_values, yend = freqs)) +
  geom_segment(aes(xend = unique_values, y = 0), color = "blue") +
  scale_y_continuous(limits = c(0, max(df$freqs))) +
  ggtitle("LMPI statistic for k = 3") +
  stat_function(fun = dnorm, args = list(mean = theta, sd = sqrt(variance)), color = "red") +
  theme_minimal()


