

load("~/nout/trials/RealData/PowerStudy/New&TidyNoSimuFunction/alpha0.2/CreditCardOnly0.2/resCredit0")
resCredit0$mean.powerGlobalNull



load("~/nout/trials/RealData/PowerStudy/New&TidyNoSimuFunction/alpha0.2/CreditCardOnly0.2/resCredit10")

boxplot(resCredit10$discoveries, main="CreditCard | Distribution of the number of discoveries")
points(x=1:5, y=resCredit10$mean.discoveries, pch=19, col="red")
resCredit10$mean.discoveries





load("~/nout/trials/RealData/PowerStudy/New&TidyNoSimuFunction/alpha0.2/DigitsOnly0.2/resDigits10")

boxplot(resDigits10$discoveries, main="Digits | Distribution of the number of discoveries")
points(x=1:5, y=resDigits10$mean.discoveries, pch=19, col="red")
resDigits10$mean.discoveries

load("~/nout/trials/RealData/PowerStudy/New&TidyNoSimuFunction/alpha0.2/DigitsOnly0.2/resDigits0")
resDigits0$mean.powerGlobalNull





load("~/nout/trials/RealData/PowerStudy/New!/alpha0.2/ShuttleOnly0.2/resShuttle10")

boxplot(resShuttle10$discoveries, main="Shuttle | Distribution of the number of discoveries")
points(x=1:5, y=resShuttle10$mean.discoveries, pch=19, col="red")
resShuttle10$mean.discoveries

load("~/nout/trials/RealData/PowerStudy/New!/alpha0.2/ShuttleOnly0.2/resShuttle0")
resShuttle0$mean.discoveries








load("~/nout/trials/RealData/PowerStudy/New!/alpha0.2/DigitsOnly0.2/Lehmann2/resDigits_alln1_alpha02_k2")

View(resDigits_alln1_alpha02_k2)

pow_WMW = vector()
pow_BH = vector()
pow_StoBH = vector()

for(i in 1:length(resDigits_alln1_alpha02_k2)){
  pow_WMW[i] = resDigits_alln1_alpha02_k2[[i]]$mean.powerGlobalNull[5]
  pow_BH[i] = resDigits_alln1_alpha02_k2[[i]]$mean.powerGlobalNull[1]
  pow_StoBH[i] = resDigits_alln1_alpha02_k2[[i]]$mean.powerGlobalNull[2]

}





load("~/nout/trials/RealData/PowerStudy/New!/alpha0.2/DigitsOnly0.2/Lehmann2/resDigits_CORRECT_alln1_alpha02_k2")


pow_WMW = vector()
pow_BH = vector()
pow_StoBH = vector()
pow_Simes = vector()
pow_StoSimes = vector()

for(i in 1:length(resDigits_alln1_alpha02_k4)){
  pow_WMW[i] = resDigits_alln1_alpha02_k4[[i]]$mean.powerGlobalNull[5]
  pow_BH[i] = resDigits_alln1_alpha02_k4[[i]]$mean.powerGlobalNull[1]
  pow_StoBH[i] = resDigits_alln1_alpha02_k4[[i]]$mean.powerGlobalNull[2]
  pow_StoSimes[i] = resDigits_alln1_alpha02_k4[[i]]$mean.powerGlobalNull[4]
  pow_Simes[i] = resDigits_alln1_alpha02_k4[[i]]$mean.powerGlobalNull[3]
}


n=50
prop.out = seq(0, 1, by=0.02)
n1_vec = round(prop.out*n)

win.graph()
plot(x = n1_vec, y = pow_WMW, col = "blue", ylab = "power",
     xlab = "n1", ylim=c(0,1), type = "b", lty = 2, pch=19,
     main = "Mean power")
points(x = n1_vec, y = pow_BH, col = "#808000", type = "b", lty = 2, pch=19)
points(x = n1_vec, y = pow_StoBH, col = "#BDB76B", type = "b", lty = 2, pch=19)
points(x = n1_vec, y = pow_Simes, col = "#DC143C", type = "b", lty = 2, pch=19)
points(x = n1_vec, y = pow_StoSimes, col = "#FFA07A", type = "b", lty = 2, pch=19)
legend("topleft", pch = 19, col = c("blue", "#808000","#BDB76B", "#DC143C", "#FFA07A"),
       legend =c("WMW", "BH", "StoreyBH", "Simes", "StoreySimes"))







pow_WMW = vector()
pow_BH = vector()
pow_StoBH = vector()
pow_Simes = vector()
pow_StoSimes = vector()

for(i in 1:length(results)){
  pow_WMW[i] = results[[i]]$mean.powerGlobalNull[5]
  pow_BH[i] = results[[i]]$mean.powerGlobalNull[1]
  pow_StoBH[i] = results[[i]]$mean.powerGlobalNull[2]
  pow_StoSimes[i] = results[[i]]$mean.powerGlobalNull[4]
  pow_Simes[i] = results[[i]]$mean.powerGlobalNull[3]
}

prop.out = seq(0, 1, by=0.02)
n1_vec = round(prop.out*n)

plot(x = n1_vec, y = pow_WMW, col = "blue", ylab = "power",
     xlab = "n1", ylim=c(0,1), type = "b", lty = 2, pch=19,
     main = "Mean power")
points(x = n1_vec, y = pow_BH, col = "#808000", type = "b", lty = 2, pch=19)
points(x = n1_vec, y = pow_StoBH, col = "#BDB76B", type = "b", lty = 2, pch=19)
points(x = n1_vec, y = pow_Simes, col = "#DC143C", type = "b", lty = 2, pch=19)
points(x = n1_vec, y = pow_StoSimes, col = "#FFA07A", type = "b", lty = 2, pch=19)
legend("topleft", pch = 19, col = c("blue", "#808000","#BDB76B", "#DC143C", "#FFA07A"),
       legend =c("WMW", "BH", "StoreyBH", "Simes", "StoreySimes"))




load("~/nout/trials/RealData/PowerStudy/New!/alpha0.2/DigitsOnly0.2/Lehmann k=4/resDigits_alpha02_k4_alln1s")
pow_WMW = vector()
pow_BH = vector()
pow_StoBH = vector()
pow_Simes = vector()
pow_StoSimes = vector()

for(i in 1:length(resDigits_alpha02_k4_alln1s)){
  pow_WMW[i] = resDigits_alpha02_k4_alln1s[[i]]$mean.powerGlobalNull[5]
  pow_BH[i] = resDigits_alpha02_k4_alln1s[[i]]$mean.powerGlobalNull[1]
  pow_StoBH[i] = resDigits_alpha02_k4_alln1s[[i]]$mean.powerGlobalNull[2]
  pow_StoSimes[i] = resDigits_alpha02_k4_alln1s[[i]]$mean.powerGlobalNull[4]
  pow_Simes[i] = resDigits_alpha02_k4_alln1s[[i]]$mean.powerGlobalNull[3]
}


n=50
prop.out = seq(0, 1, by=0.02)
n1_vec = round(prop.out*n)

win.graph()
plot(x = n1_vec, y = pow_WMW, col = "blue", ylab = "power",
     xlab = "n1", ylim=c(0,1), type = "b", lty = 2, pch=19,
     main = "Mean power")
points(x = n1_vec, y = pow_BH, col = "#808000", type = "b", lty = 2, pch=19)
points(x = n1_vec, y = pow_StoBH, col = "#BDB76B", type = "b", lty = 2, pch=19)
points(x = n1_vec, y = pow_Simes, col = "#DC143C", type = "b", lty = 2, pch=19)
points(x = n1_vec, y = pow_StoSimes, col = "#FFA07A", type = "b", lty = 2, pch=19)
legend("topleft", pch = 19, col = c("blue", "#808000","#BDB76B", "#DC143C", "#FFA07A"),
       legend =c("WMW", "BH", "StoreyBH", "Simes", "StoreySimes"))






load("~/nout/trials/RealData/PowerStudy/New!/alpha0.2/DigitsOnly0.2/Lehmann k=4/resDigits_alpha02_k4_alln1s")
pow_WMW = vector()
pow_BH = vector()
pow_StoBH = vector()
pow_Simes = vector()
pow_StoSimes = vector()

for(i in 1:length(resDigits_alpha02_k4_alln1s)){
  pow_WMW[i] = resDigits_alpha02_k4_alln1s[[i]]$mean.powerGlobalNull[5]
  pow_BH[i] = resDigits_alpha02_k4_alln1s[[i]]$mean.powerGlobalNull[1]
  pow_StoBH[i] = resDigits_alpha02_k4_alln1s[[i]]$mean.powerGlobalNull[2]
  pow_StoSimes[i] = resDigits_alpha02_k4_alln1s[[i]]$mean.powerGlobalNull[4]
  pow_Simes[i] = resDigits_alpha02_k4_alln1s[[i]]$mean.powerGlobalNull[3]
}


n=50
prop.out = seq(0, 1, by=0.02)
n1_vec = round(prop.out*n)

win.graph()
plot(x = n1_vec, y = pow_WMW, col = "blue", ylab = "power",
     xlab = "n1", ylim=c(0,1), type = "b", lty = 2, pch=19,
     main = "Mean power")
points(x = n1_vec, y = pow_BH, col = "#808000", type = "b", lty = 2, pch=19)
points(x = n1_vec, y = pow_StoBH, col = "#BDB76B", type = "b", lty = 2, pch=19)
points(x = n1_vec, y = pow_Simes, col = "#DC143C", type = "b", lty = 2, pch=19)
points(x = n1_vec, y = pow_StoSimes, col = "#FFA07A", type = "b", lty = 2, pch=19)
legend("topleft", pch = 19, col = c("blue", "#808000","#BDB76B", "#DC143C", "#FFA07A"),
       legend =c("WMW", "BH", "StoreyBH", "Simes", "StoreySimes"))





load("~/nout/trials/RealData/PowerStudy/New!/alpha0.2/CreditCardOnly0.2/Lehmann k=2/resCreditCard_alpha02_k2_alln1s")
pow_WMW = vector()
pow_BH = vector()
pow_StoBH = vector()
pow_Simes = vector()
pow_StoSimes = vector()

for(i in 1:length(resDigits_alpha02_k2_alln1s)){
  pow_WMW[i] = resDigits_alpha02_k2_alln1s[[i]]$mean.powerGlobalNull[5]
  pow_BH[i] = resDigits_alpha02_k2_alln1s[[i]]$mean.powerGlobalNull[1]
  pow_StoBH[i] = resDigits_alpha02_k2_alln1s[[i]]$mean.powerGlobalNull[2]
  pow_StoSimes[i] = resDigits_alpha02_k2_alln1s[[i]]$mean.powerGlobalNull[4]
  pow_Simes[i] = resDigits_alpha02_k2_alln1s[[i]]$mean.powerGlobalNull[3]
}


n=50
prop.out = seq(0, 1, by=0.02)
n1_vec = round(prop.out*n)

win.graph()
plot(x = n1_vec, y = pow_WMW, col = "blue", ylab = "power",
     xlab = "n1", ylim=c(0,1), type = "b", lty = 2, pch=19,
     main = "Credit mean power | k=2")
points(x = n1_vec, y = pow_BH, col = "#808000", type = "b", lty = 2, pch=19)
points(x = n1_vec, y = pow_StoBH, col = "#BDB76B", type = "b", lty = 2, pch=19)
points(x = n1_vec, y = pow_Simes, col = "#DC143C", type = "b", lty = 2, pch=19)
points(x = n1_vec, y = pow_StoSimes, col = "#FFA07A", type = "b", lty = 2, pch=19)
legend("topleft", pch = 19, col = c("blue", "#808000","#BDB76B", "#DC143C", "#FFA07A"),
       legend =c("WMW", "BH", "StoreyBH", "Simes", "StoreySimes"))





