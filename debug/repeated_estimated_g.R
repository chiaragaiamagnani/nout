# ------------------------------------------------------------- #
# PLOTS OF THE POWER OF WMW, FISHER, ORACLE AND Shiraishi TESTS #
# ------------------------------------------------------------- #

rm(list = ls())


library(tidyverse)
library(latex2exp)
library(gridExtra)
library(foreach)
library(tictoc)
library(fitdistrplus)

library(nout)
source("~/nout/R/utils_g.R")
source("~/nout/R/utils_higher.R")
source("~/nout/R/utils_global_functions.R")
source("~/nout/debug/utils_g_simu_repeated.R")
source("~/nout/R/utils_LehmannAlt.R")


#### Lehmann's alternative corresponding to k=2,3,10 (k>=2) ####


B=100; B_MC=500 # number of repetitions
sizes = c(150) # m=n
thetas = seq(from=0, to=1, by=0.05) # proportion of outliers in the test set
alpha = 0.1 # significance level


# k values and corresponding distributions (null and alternative)

ks = c(2,3,10)
rg_null = function(x){ runif(x) }
rg_alt_list = list(function(x){ rg1(x, rg_null=rg_null) },
                   function(x){ rg2(x, rg_null=rg_null) },
                   function(x){ rg9(x, rg_null=rg_null) })
g_alt_list = list(function(x){ g1(x) },
                  function(x){ g2(x) },
                  function(x){ g9(x) })


# Run simulation

tic()
res = lapply(1:length(rg_alt_list),
             function(ll){
               cat("k=",ks[ll],"\n")
               res3 = lapply(1:length(sizes), function(s){
                 cat("n=",sizes[s],"\n")
                 res2 = sapply(1:length(thetas),
                               function(j){
                                 cat("theta=",thetas[j],"\n")
                                 res1 = sim_pow_rep(B=B, B_MC=B_MC, m=sizes[s], n=sizes[s], theta=thetas[j],
                                                rg_null=rg_null, rg=rg_alt_list[[ll]], g=g_alt_list[[ll]],
                                                alpha=alpha)
                                 out1 = tibble(k=ks[ll],theta=thetas[j],res1)
                               })

               })
               out2 = res3[[1]]
               # To be commented if length(sizes) = 1
               # for(i in 2:length(res3)){
               #   out2 = cbind(out2,res3[[i]])
               # }
               return(out2)

             })
toc()

out = res[[1]]
for(i in 2:length(res)){
  out = cbind(out,res[[i]])
}
out = t(out)
c.names = colnames(out)
out = as.data.frame(sapply(1:ncol(out), function(j) unlist(out[,j])))
colnames(out) = c.names

View(out)

# pow_LehmannAlt_B1000_B_MC1000_mixprop2016_TFBoth = out
# save(pow_LehmannAlt_B1000_B_MC1000_mixprop2016_TFBoth,
#     file = "~/nout/mynotes/PowerComparison/pow_LehmannAlt_B1000_B_MC1000_mixprop2016_TFBoth")

pow_LehmannAlt_B100_B_MC500 = out
save(pow_LehmannAlt_B100_B_MC500,
     file = "~/nout/mynotes/PowerComparison/pow_LehmannAlt_B100_B_MC500")


# Plots

# load("~/nout/mynotes/PowerComparison/paper_pow_v2_LehmannAlt_B500_B_MC100")
# out = paper_pow_v2_LehmannAlt_B500_B_MC100


method.values <- c("WMW","Fisher","oracle", "KDE", "mono")
method.labels <- c("WMW","Fisher","oracle", "KDE", "mono")
my_colors <- c("#377EB8","#37b8b2", "#F3C608","#91099e", "#A87AAD") # "#984EA3"  EA5B17

df = as_tibble(out) %>%
  as.data.frame %>%
  mutate(n = factor(n)) %>%
  pivot_longer(c("power_oracle", "power_WMW", "power_Fisher", "power_KDE", "power_mono",
                 "std.error_oracle", "std.error_WMW", "std.error_Fisher", "std.error_KDE", "std.error_mono"),
               names_to = c(".value", "Method"), names_sep = "_") %>%
  mutate(power = as.double(power))  %>%
  mutate(Method = factor(Method, method.values, method.labels))#%>%
#pivot_longer(c("power", "std.error"), names_to="Key", values_to="Value")


pp <- df %>%
  ggplot(aes(x=theta, y=power, color=Method, alpha=Method, shape=Method)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = power - std.error,
                    ymax = power + std.error),
                width = 0.02,
                position = position_dodge(width = 0.013),
                alpha = 0.9) +
  facet_grid(n ~ k,
             labeller = labeller(n = label_both, k = label_both),
             scales = "free")+
  # facet_wrap(~n+k, labeller=labeller(n=label_both,k = label_both), scales="free")+
  scale_color_manual(values=my_colors) +
  #scale_color_brewer(palette = "Set1") +
  scale_alpha_manual(values = c(1,1,1,1.5,1.5,1.5)) +
  xlab(TeX("\\theta")) +
  ylab("") +
  xlim(-0.01,max(df$theta)+0.01) +
  ylim(-0.1,1) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  guides(color=guide_legend(nrow=1, byrow=TRUE))+
  ggtitle("Lehmann's alternatives")

pp




#### Location shift with normal null distribution ####

B=100; B_MC=500 # number of repetitions
sizes = c(10,100,1000) # m=n
thetas = seq(from=0, to=0.2, by=0.05) # proportion of outliers in the test set
alpha = 0.1 # significa level


# shift values

shift = c(-1,1)

rg_null = function(x){ rnorm(x) }
rg_alt_list = list(function(x){ rnorm(n=x, mean=shift[1]) },
                   function(x){ rnorm(n=x, mean=shift[2]) })

g_alt_loc_shift = function(x,mu)   dnorm( qnorm(x) - mu )*(1/dnorm(qnorm(x)))

g_alt_list = list(function(x){ g_alt_loc_shift(x=x, mu=shift[1]) },
                  function(x){ g_alt_loc_shift(x=x, mu=shift[2]) })


# Add a for loop so that for each alternative distribution g is estimated

X = rg_null(x=4000)
X1 = X[1:2000]
X2 = X[2001:4000]

g_hat = list()
g_hat_monotone = list()
constr = c("decreasing", "increasing")

tic()
for(i in 1:length(g_alt_list)){
  Y = rg_alt_list[[i]](x=2000)
  g_hat_monotone[[i]] = estimate_g(X1=X1, X2=X2, Y=Y, constraint=constr[i], ker="uniform")
}
toc()


tic()
for(i in 1:length(g_alt_list)){
  Y = rg_alt_list[[i]](x=2000)
  g_hat[[i]] = estimate_g(X1=X1, X2=X2, Y=Y, ker="uniform")
}
toc()



# Run simulation

tic()
res = lapply(1:length(rg_alt_list),
             function(ll){
               cat("shift=",shift[ll],"\n")
               res3 = lapply(1:length(sizes), function(s){
                 cat("n=",sizes[s],"\n")
                 res2 = sapply(1:length(thetas),
                               function(j){
                                 cat("theta=",thetas[j],"\n")
                                 res1 = sim_pow_rep(B=B, B_MC=B_MC, m=sizes[s], n=sizes[s], theta=thetas[j],
                                                rg_null=rg_null, rg=rg_alt_list[[ll]], g=g_alt_list[[ll]],
                                                g_hat=g_hat[[ll]], g_hat_monotone=g_hat_monotone[[ll]],
                                                alpha=alpha)
                                 out1 = tibble(shift=shift[ll],theta=thetas[j],res1)
                               })

               })
               out2 = res3[[1]]
               for(i in 2:length(res3)){
                 out2 = cbind(out2,res3[[i]])
               }
               return(out2)

             })
toc()


out = res[[1]]
for(i in 2:length(res)){
  out = cbind(out,res[[i]])
}
out = t(out)
c.names = colnames(out)
out = as.data.frame(sapply(1:ncol(out), function(j) unlist(out[,j])))
colnames(out) = c.names

View(out)


# pow_LocationShift_B1000_B_MC1000_mixpropPS2016_TFBoth = out
# save(pow_LocationShift_B1000_B_MC1000_mixpropPS2016_TFBoth,
#      file = "~/nout/mynotes/PowerComparison/pow_LocationShift_B1000_B_MC1000_mixpropPS2016_TFBoth")

pow_LocationShift_B100_B_MC500 = out
save(pow_LocationShift_B100_B_MC500,
     file = "~/nout/mynotes/PowerComparison/pow_LocationShift_B100_B_MC500")


# Plots

# load("~/nout/mynotes/PowerComparison/pow_LocationShift_B1000_B_MC1000_mixpropPS2016_TFBoth")
# out = pow_LocationShift_B1000_B_MC1000_mixpropPS2016_TFBoth

method.values <- c("WMW","Fisher","oracle", "estG", "estGmono")
method.labels <- c("WMW","Fisher","oracle", "estG", "estGmono")
my_colors <- c("#377EB8","#37b8b2", "#F3C608","#91099e", "#A87AAD")

df = as_tibble(out) %>%
  as.data.frame %>%
  mutate(n = factor(n)) %>%
  pivot_longer(c("power_oracle", "power_WMW", "power_Fisher", "power_estG", "power_estGmono",
                 "std.error_oracle", "std.error_WMW", "std.error_Fisher", "std.error_estG", "std.error_estGmono"),
               names_to = c(".value", "Method"), names_sep = "_") %>%
  mutate(power = as.double(power))  %>%
  mutate(Method = factor(Method, method.values, method.labels))#%>%
#pivot_longer(c("power", "std.error"), names_to="Key", values_to="Value")


pp <- df %>%
  ggplot(aes(x=theta, y=power, color=Method, alpha=Method, shape=Method)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = power - std.error,
                    ymax = power + std.error),
                width = 0.02,
                position = position_dodge(width = 0.013),
                alpha = 0.9) +
  facet_grid(shift ~ n,
             labeller = labeller(n = label_both, shift = label_both),
             scales = "free")+
  # facet_wrap(~n+k, labeller=labeller(n=label_both,k = label_both), scales="free")+
  scale_color_manual(values = my_colors) +
  scale_alpha_manual(values = c(1,0.5,0.5,0.5,1,1)) +
  xlab(TeX("\\theta")) +
  ylab("") +
  xlim(-0.01,max(df$theta)+0.01) +
  ylim(-0.1,1) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  guides(color=guide_legend(nrow=1, byrow=TRUE))+
  ggtitle("Location shift")

pp






#### Concentration difference with normal null distribution ####

B=100; B_MC=500 # number of repetitions
sizes = c(10,100,1000) # m=n
thetas = seq(from=0, to=0.2, by=0.05) # proportion of outliers in the test set
alpha = 0.1 # significa level


# sd values
sd = c(0.2, 0.9, 1.8)
sd = c(0.2, 1.8)

rg_null = function(x){ rnorm(x) }
rg_alt_list = list(function(x){ rnorm(x, sd=0.2) },
                   function(x){ rnorm(x, sd=1.8) })

g_alt_scale = function(x,sigma)   dnorm( qnorm(x)/sigma )*(1/(sigma*dnorm(qnorm(x))))

g_alt_list = list(function(x){ g_alt_scale(x=x, sigma=sd[1]) },
                  function(x){ g_alt_scale(x=x, sigma=sd[2]) })



# Add a for loop so that for each alternative distribution g is estimated

X = rg_null(x=4000)
X1 = X[1:2000]
X2 = X[2001:4000]

g_hat = list()
g_hat_monotone = list()

tic()
for(i in 1:length(g_alt_list)){
  Y = rg_alt_list[[i]](x=2000)
  g_hat_monotone[[i]] = estimate_g(X1=X1, X2=X2, Y=Y, constraint="increasing", ker="uniform")
}
toc()


tic()
for(i in 1:length(g_alt_list)){
  Y = rg_alt_list[[i]](x=2000)
  g_hat[[i]] = estimate_g(X1=X1, X2=X2, Y=Y, ker="uniform")
}
toc()

# Run simulation

tic()
res = lapply(1:length(rg_alt_list),
             function(ll){
               cat("sd=",sd[ll],"\n")
               res3 = lapply(1:length(sizes), function(s){
                 cat("n=",sizes[s],"\n")
                 res2 = sapply(1:length(thetas),
                               function(j){
                                 cat("theta=",thetas[j],"\n")
                                 res1 = sim_pow_rep(B=B, B_MC=B_MC, m=sizes[s], n=sizes[s], theta=thetas[j],
                                                rg_null=rg_null, rg=rg_alt_list[[ll]], g=g_alt_list[[ll]],
                                                g_hat=g_hat[[ll]], g_hat_monotone=g_hat_monotone[[ll]],
                                                alpha=alpha)
                                 out1 = tibble(sd=sd[ll],theta=thetas[j],res1)
                               })

               })
               out2 = res3[[1]]
               for(i in 2:length(res3)){
                 out2 = cbind(out2,res3[[i]])
               }
               return(out2)

             })
toc()


out = res[[1]]
for(i in 2:length(res)){
  out = cbind(out,res[[i]])
}
out = t(out)
c.names = colnames(out)
out = as.data.frame(sapply(1:ncol(out), function(j) unlist(out[,j])))
colnames(out) = c.names

View(out)

# pow_sdNormAlt_B1000_B_MC1000 = out
# save(pow_sdNormAlt_B1000_B_MC1000,
#      file = "~/nout/mynotes/PowerComparison/pow_sdNormAlt_B1000_B_MC1000")

paper_pow_v2_sdNormAlt_B500_B_MC1000 = out
save(paper_pow_v2_sdNormAlt_B500_B_MC1000,
     file = "~/nout/mynotes/PowerComparison/paper_pow_v2_sdNormAlt_B500_B_MC1000")


# Plots

# load("~/nout/mynotes/PowerComparison/pow_sdNormAlt_B1000_B_MC1000")
# load("~/nout/mynotes/PowerComparison/paper_pow_sdNormAlt_B1000_B_MC1000_mixprop2016")
# out = paper_pow_sdNormAlt_B1000_B_MC1000_mixprop2016

method.values <- c("WMW","Fisher","oracle", "estG", "estGmono")
method.labels <- c("WMW","Fisher","oracle", "estG", "estGmono")
my_colors <- c("#377EB8","#37b8b2", "#F3C608","#750880", "#A87AAD")

df = as_tibble(out) %>%
  as.data.frame %>%
  mutate(n = factor(n)) %>%
  pivot_longer(c("power_oracle", "power_WMW", "power_Fisher", "power_estG", "power_estGmono",
                 "std.error_oracle", "std.error_WMW", "std.error_Fisher", "std.error_estG", "std.error_estGmono"),
               names_to = c(".value", "Method"), names_sep = "_") %>%
  mutate(power = as.double(power))  %>%
  mutate(Method = factor(Method, method.values, method.labels))


pp <- df %>%
  ggplot(aes(x=theta, y=power, color=Method, alpha=Method, shape=Method)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = power - std.error,
                    ymax = power + std.error),
                width = 0.02,
                position = position_dodge(width = 0.013),
                alpha = 0.9) +
  facet_grid(sd ~ n,
             labeller = labeller(n = label_both, shift = label_both),
             scales = "free")+
  # facet_wrap(~n+k, labeller=labeller(n=label_both,k = label_both), scales="free")+
  scale_color_manual(values = my_colors) +
  # scale_color_brewer(palette = "Set1") +
  scale_alpha_manual(values = c(1,0.5,0.5,0.5, 1, 1)) +
  xlab(TeX("\\theta")) +
  ylab("") +
  xlim(-0.01,max(df$theta)+0.01) +
  ylim(-0.1,1) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  guides(color=guide_legend(nrow=1, byrow=TRUE))+
  ggtitle("Different scale")

pp






#### Beta alternatives ####

# different parameters values and corresponding distributions (null and alternative)

B=100; B_MC=500 # number of repetitions
sizes = c(10,100,1000) # m=n
thetas = seq(from=0, to=0.2, by=0.05) # proportion of outliers in the test set
alpha = 0.1 # significa level

params = as.data.frame(rbind(c(0.25,0.25),c(25,25),c(10,1),c(500,300)))
params = as.data.frame(rbind(c(0.25,0.25),c(25,25),c(500,300)))

colnames(params) = c("a","b")

rg_null = function(x){ runif(x) }

# rg_alt_list = list(function(x){ rbeta(n=x, shape1 = params$a[1], shape2 = params$b[1]) }, # overdispersion
#                    function(x){ rbeta(n=x, shape1 = params$a[2], shape2 = params$b[2]) }, # underdispersion
#                    function(x){ rbeta(n=x, shape1 = params$a[3], shape2 = params$b[3]) }, # Lehmann
#                    function(x){ rbeta(n=x, shape1 = params$a[4], shape2 = params$b[4]) }) # weak concentrated effects
#
# g_alt_list = list(function(x){ dbeta(x, shape1 = params$a[1], shape2 = params$b[1]) },
#                   function(x){ dbeta(x, shape1 = params$a[2], shape2 = params$b[2]) },
#                   function(x){ dbeta(x, shape1 = params$a[3], shape2 = params$b[3]) },
#                   function(x){ dbeta(x, shape1 = params$a[4], shape2 = params$b[4]) })

rg_alt_list = list(function(x){ rbeta(n=x, shape1 = params$a[1], shape2 = params$b[1]) }, # overdispersion
                   function(x){ rbeta(n=x, shape1 = params$a[2], shape2 = params$b[2]) }, # underdispersion
                   function(x){ rbeta(n=x, shape1 = params$a[3], shape2 = params$b[3]) }) # weak concentrated effects

g_alt_list = list(function(x){ dbeta(x, shape1 = params$a[1], shape2 = params$b[1]) },
                  function(x){ dbeta(x, shape1 = params$a[2], shape2 = params$b[2]) },
                  function(x){ dbeta(x, shape1 = params$a[3], shape2 = params$b[3]) })

# rg_alt_list = list(function(x){ rbeta(n=x, shape1 = params$a[1], shape2 = params$b[1]) }) # weak concentrated effects
# g_alt_list = list(function(x){ dbeta(x, shape1 = params$a[1], shape2 = params$b[1]) })

# Add a for loop so that for each alternative distribution g is estimated

X = rg_null(x=4000)
X1 = X[1:2000]
X2 = X[2001:4000]

g_hat = list()
g_hat_monotone = list()

tic()
for(i in 1:length(g_alt_list)){
  Y = rg_alt_list[[i]](x=2000)
  g_hat_monotone[[i]] = estimate_g(X1=X1, X2=X2, Y=Y, constraint="increasing", ker="uniform")
}
toc()


tic()
for(i in 1:length(g_alt_list)){
  Y = rg_alt_list[[i]](x=2000)
  g_hat[[i]] = estimate_g(X1=X1, X2=X2, Y=Y, ker="uniform")
}
toc()


# Run experiment
tic()
res = lapply(1:length(rg_alt_list),
             function(ll){
               cat("a=",params$a[ll],"\n")
               res3 = lapply(1:length(sizes), function(s){
                 cat("n=",sizes[s],"\n")
                 res2 = sapply(1:length(thetas),
                               function(j){
                                 cat("theta=",thetas[j],"\n")
                                 res1 = sim_pow_rep(B=B, B_MC=B_MC, m=sizes[s], n=sizes[s], theta=thetas[j],
                                                rg_null=rg_null, rg=rg_alt_list[[ll]], g=g_alt_list[[ll]],
                                                g_hat=g_hat[[ll]], g_hat_monotone=g_hat_monotone[[ll]],
                                                alpha=alpha)
                                 out1 = tibble(a=params$a[ll],b=params$b[ll],theta=thetas[j],res1)
                               })

               })
               out2 = res3[[1]]
               for(i in 2:length(res3)){
                 out2 = cbind(out2,res3[[i]])
               }
               return(out2)

             })
toc()


out = res[[1]]
for(i in 2:length(res)){
  out = cbind(out,res[[i]])
}
out = t(out)
c.names = colnames(out)
out = as.data.frame(sapply(1:ncol(out), function(j) unlist(out[,j])))
colnames(out) = c.names

View(out)

# pow_BetaAlt_B1000_B_MC1000 = out
# save(pow_BetaAlt_B1000_B_MC1000,
#     file = "~/nout/mynotes/PowerComparison/pow_BetaAlt_B1000_B_MC1000")

pow_BetaAlt_B100_B_MC500 = out
save(pow_BetaAlt_B100_B_MC500,
     file = "~/nout/mynotes/PowerComparison/pow_BetaAlt_B100_B_MC500")




# Plots

# load("~/nout/mynotes/PowerComparison/pow_BetaAlt_B1000_B_MC1000")
# load("~/nout/mynotes/PowerComparison/pow_BetaAlt_B1000_B_MC1000_mixprop2016_TFBoth")
# out = pow_BetaAlt_B1000_B_MC1000_mixprop2016_TFBoth

method.values <- c("WMW","Fisher","oracle", "estG", "estGmono")
method.labels <- c("WMW","Fisher","oracle", "estG", "estGmono")
my_colors <- c("#377EB8","#37b8b2", "#F3C608","#750880", "#A87AAD") #"#91099e"

df = as_tibble(out) %>%
  as.data.frame %>%
  mutate(n = factor(n)) %>%
  pivot_longer(c("power_oracle", "power_WMW", "power_Fisher", "power_estG", "power_estGmono",
                 "std.error_oracle", "std.error_WMW", "std.error_Fisher", "std.error_estG", "std.error_estGmono"),
               names_to = c(".value", "Method"), names_sep = "_") %>%
  mutate(power = as.double(power))  %>%
  mutate(Method = factor(Method, method.values, method.labels))


pp <- df %>%
  ggplot(aes(x=theta, y=power, color=Method, alpha=Method, shape=Method)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = power - std.error,
                    ymax = power + std.error),
                width = 0.02,
                position = position_dodge(width = 0.013),
                alpha = 0.9) +
  facet_grid(a+b ~ n,
             labeller = labeller(n = label_both, a = label_both, b = label_both),
             scales = "free")+
  # facet_wrap(~n+k, labeller=labeller(n=label_both,k = label_both), scales="free")+
  scale_color_manual(values=my_colors) +
  scale_alpha_manual(values = c(1,0.8,0.7,0.6,0.5,0.4)) +
  xlab(TeX("\\theta")) +
  ylab("") +
  xlim(-0.01,max(df$theta)+0.01) +
  ylim(-0.1,1) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  guides(color=guide_legend(nrow=1, byrow=TRUE))+
  ggtitle("Beta Alternatives")

pp





