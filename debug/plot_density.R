
# Plot (estimated) densities

df <- as.data.frame(res1$dens)
colnames(df) = c("u", "g_true", "g_KDE", "g_mono")

method.values = c("true", "KDE", "mono")
method.labels = c("true", "KDE", "mono")
my_colors = c("black", "red", "blue")
df = as_tibble(df) %>%
  as.data.frame %>%
  pivot_longer(c("g_true", "g_KDE", "g_mono"),
               names_to = c(".value", "Method"), names_sep = "_") %>%
  mutate(Method = factor(Method, method.values, method.labels))

pp1 = df %>%
  ggplot(aes(x=u, y=g, color=Method, alpha=Method, shape=Method)) +
  geom_line(size=1) +
  scale_color_manual(values=my_colors) +
  scale_alpha_manual(values = c(1,0.8,0.7)) +
  xlab("u") +
  ylab("g") +
  theme_bw() +
  theme(legend.position="right",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  guides(color=guide_legend(nrow=3, byrow=TRUE))






u = runif(5000)
g.true = g_alt_list[[1]](u)
g.KDE = g_KDE(u)
g.mono = g_mono(u)
g.KDE2 = g_KDE2(u)
g.mono2 = g_mono2(u)

df <- as.data.frame(cbind(u, g.true, g.KDE, g.mono, g.KDE2, g.mono2))
colnames(df) = c("u", "g_true", "g_KDE", "g_mono", "g_KDE2", "g_mono2")

method.values = c("true", "KDE", "mono", "KDE2", "mono2")
method.labels = c("true", "KDE", "mono", "KDE2", "mono2")
my_colors = c("black", "red", "blue", "pink", "lightblue")
df = as_tibble(df) %>%
  as.data.frame %>%
  pivot_longer(c("g_true", "g_KDE", "g_mono", "g_KDE2", "g_mono2"),
               names_to = c(".value", "Method"), names_sep = "_") %>%
  mutate(Method = factor(Method, method.values, method.labels))

pp1 = df %>%
  ggplot(aes(x=u, y=g, color=Method, alpha=Method, shape=Method)) +
  geom_line(size=1) +
  scale_color_manual(values=my_colors) +
  scale_alpha_manual(values = c(1,0.8,0.7, 0.6, 0.5)) +
  xlab("u") +
  ylab("g") +
  theme_bw() +
  theme(legend.position="right",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  guides(color=guide_legend(nrow=3, byrow=TRUE))

pp1













sim_Tk = function(m,n,k){
  X = runif(m)
  Y = runif(n)
  Tk = sum(stat.Tk(Z=c(X,Y), m=m, k=k))
  return(Tk)
}

m=100
n=100
k=1
stats=replicate(5000, sim_Tk(m=m, n=n, k=k))

hist(stats, freq = F)
mom = asymptotic.moments.Tk(m=m, n=n, k=k)
curve(dnorm(x, mean=mom$mean.Tk, sd = sqrt(mom$variance.Tk)), add=T, col="red")




