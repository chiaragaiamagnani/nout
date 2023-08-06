load("~/nout/trials/RealData/PowerStudy/FinalSimu/Digits/resDigits0.1")

disc_Sim = vector()
disc_StoSimes = vector()
disc_WMW = vector()
disc_WMW.cpp = vector()

for(j in 1:length(n1s)){
  disc_Sim[j] = resDigits0.1[["compact.results"]][[j]]$mean.n.disc[1]
  disc_StoSimes[j] = resDigits0.1[["compact.results"]][[j]]$mean.n.disc[3]
  disc_WMW[j] = resDigits0.1[["compact.results"]][[j]]$mean.n.disc[4]
  disc_WMW.cpp[j] = resDigits0.1[["compact.results"]][[j]]$mean.n.disc[5]
}

df <- data.frame(
  x = n1s,
  Simes_CT = disc_Sim,
  StoreySimes_CT = disc_StoSimes,
  WMW.cpp_CT = disc_WMW.cpp,
  WMW_CT = disc_WMW
)
df_long <- tidyr::pivot_longer(df, cols = -x, names_to = "group", values_to = "y")

ggplot(df_long, aes(x = x, y = y, color = group)) +
  geom_line() +
  geom_point()+
  scale_color_manual(values = c("#DC143C", "#FFA07A", "blue", 5)) +
  labs(x = "n1", y = "d", title = "Mean of the number of discoveries on B replications") +
  theme_minimal() +
  theme(legend.title = element_blank())

