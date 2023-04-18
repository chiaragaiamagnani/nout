
m=10;n=20
crits = critWMW(m=m, n=n)$crit.vals
plot(x=1:m, y=crits, xlab="m-k+1", main = "Critical values of WMW statistic for n=20 and k=1,...,m")



plot(x=1:(m*n),y=dwilcox(1:(m*n),m=m,n=n), ylim=c(0,0.05),
     main = "Density function of WMW statistic for n=20")
lines(dwilcox(1:(m*n),m=m-1,n=n), lwd = 2, col = "violet")
lines(dwilcox(1:(m*n),m=m-2,n=n), lwd = 2, col = "blue")
lines(dwilcox(1:(m*n),m=m-3,n=n), lwd = 2, col = "yellow")
lines(dwilcox(1:(m*n),m=m-4,n=n), lwd = 2, col = "green")
lines(dwilcox(1:(m*n),m=m-5,n=n), lwd = 2, col = "orange")
lines(dwilcox(1:(m*n),m=m-6,n=n), lwd = 2, col = "red")
lines(dwilcox(1:(m*n),m=m-7,n=n), lwd = 2, col = "lightblue")
lines(dwilcox(1:(m*n),m=m-8,n=n), lwd = 2, col = "darkgreen")
lines(dwilcox(1:(m*n),m=m-9,n=n), lwd = 2, col = "purple")
legend("topright", legend = c("m=10", "m=9", "m=8", "m=7",
                              "m=6", "m=5", "m=4", "m=3", "m=2", "m=1"),
       col = c("black", "violet", "blue", "yellow", "green", "orange", "red",
               "lightblue", "darkgreen", "purple"), lty=1)
