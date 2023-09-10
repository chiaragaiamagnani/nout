
# Comportamento asintotico Wilcoxon Mann-Whitney
m = 100
n = 100
B = 10^4

U = vector()

for(b in 1:B){

  X = runif(n=m, min=0, max=2)
  Y = runif(n=n, min=0, max=2)

  Uj = sapply(1:n, function(j) sum(X<Y[j]))
  U[b] = sum(Uj)

}

var(U)
plot(sapply(1:B, function(b) var(U[1:b])))
plot(sapply(1:B, function(b) mean(U[1:b])))
theta = m*n/2
varr = (m*n*(m+n))/12



# Comportamento asintotico statistica test per k=3
m = 400
n = 400
B = 10^4

U = vector()

for(b in 1:B){

  X = runif(n=m, min=0, max=2)
  Y = runif(n=n, min=0, max=2)

  Uj = sapply(1:n, function(j) sum(X<Y[j]))
  Uj2 = Uj^2
  U[b] = sum(Uj2+Uj)

}
mean(U)
# var(U)
# plot(sapply(1:B, function(b) var(U[1:b])))
plot(sapply(1:B, function(b) mean(U[1:b])))
(theta = n*(m*0.5+m^2/3))
# varr







# Verifico correttezza della riscrittura della funzione kernel
m = 100
n = 100
B = 30

Uv1 = vector()
Uv2 = vector()

for(b in 1:B){

  X = runif(n=m, min=0, max=2)
  Y = runif(n=n, min=0, max=2)

  Uj = sapply(1:n, function(j) sum(X<Y[j]))
  Uj2 = Uj^2
  Uv1[b] = sum(Uj2+Uj)

  cond = matrix(nrow=n, ncol=m)
  for(j in 1:n){
    for(i in 1:m){
      if(X[i]<Y[j]){
        for(h in 1:m){
          cond[j,i] = sum(X<Y[j])
        }
      }
      else{
        cond[j,i] = F
      }
    }
  }
  condd = apply(cond, MARGIN = 1, FUN = sum)
  Uv2[b] = sum(condd+Uj)

}

Uv2==Uv1




# Verifico correttezza valore atteso della funzione kernel h
ciclefun = function(indices, index, cond){
  res = sapply(indices, function(j) min(cond[j],cond[index]))
  return(res)
}




hfun = function(x, y){
  cond = x<y
  aus = vector()
  nr = length(x)

  condoublematrix = sapply(1:nr, function(i) ciclefun(indices=1:nr, index=i, cond=cond))

  h = sum(cond)+sum(apply(condoublematrix, MARGIN = 1, FUN=sum))
  #h1 = sum(cond)+sum(sapply(1:nr, function(i) sum(condoublematrix[i,])))

  return(h)
}



m = 10
n = 10
B = 10^7
set.seed(123)
kernelvalues = vector()
X1 = runif(n=1, min=0, max=2)

for(b in 1:B){

  X = c(X1, runif(n=(m-1), min=0, max=2))
  Y = runif(n=1, min=0, max=2)

  kernelvalues[b] = hfun(X, Y)

}
mean(kernelvalues)

(2*m^2+5*m+5)/6-2*punif(X1, min = 0, max = 2)-(m-1)*(punif(X1, min = 0, max = 2))^2






# Verifico esattezza della forma chiusa del valore atteso di 1{x1<Y, X<Y}, dove Y,X random e X1 fissato

B=10^6

X1 = runif(n=1, min=0, max=2)
res = vector()

for(b in 1:B){

  X = runif(n=1, min=0, max=2)
  Y = runif(n=1, min=0, max=2)

  res[b] = ifelse((X1<Y & X<Y), 1, 0)

}
mean(res)
0.5*(1-punif(X1, min = 0, max = 2))^2+(1-punif(X1, min = 0, max = 2))*punif(X1, min = 0, max = 2)




# Verifico esattezza della forma chiusa del valore atteso di 1{x1<Y, X2<Y}, dove Y,X1,X2 random
B=10^6

res = vector()

for(b in 1:B){

  X1 = runif(n=1, min=0, max=2)
  X2 = runif(n=1, min=0, max=2)
  Y = runif(n=1, min=0, max=2)

  res[b] = ifelse((X1<Y & X2<Y), 1, 0)

}
mean(res)
1/3



# Verifico esattezza della forma chiusa del valore atteso di 1{x<Y}, dove Y,X random
B=10^6

res = vector()

for(b in 1:B){

  X = runif(n=1, min=0, max=2)
  Y = runif(n=1, min=0, max=2)

  res[b] = ifelse(X<Y, 1, 0)

}
mean(res)
1/2




# Verifico esattezza della forma chiusa di P(x1<Xh<Xr)

B=10^6

res = vector()
x1 = runif(n=1, min=0, max=2)
for(b in 1:B){

  Xh = runif(n=1, min=0, max=2)
  Xr = runif(n=1, min=0, max=2)

  res[b] = ifelse(Xh<Xr & x1<Xh, 1, 0)

}
mean(res)
0.5*(1-punif(x1, min = 0, max = 2))^2





# Verifico esattezza della forma chiusa di P(Xh<x1<Xr)

B=10^6

res = vector()
x1 = runif(n=1, min=0, max=2)
for(b in 1:B){

  Xh = runif(n=1, min=0, max=2)
  Xr = runif(n=1, min=0, max=2)

  res[b] = ifelse(Xh<x1 & x1<Xr, 1, 0)

}
mean(res)
(1-punif(x1, min = 0, max = 2))*(punif(x1, min = 0, max = 2))






# Verifico esattezza della forma chiusa di P(x1>Xh>Xr)

B=10^6

res = vector()
x1 = runif(n=1, min=0, max=2)
for(b in 1:B){

  Xh = runif(n=1, min=0, max=2)
  Xr = runif(n=1, min=0, max=2)

  res[b] = ifelse(x1>Xh & Xh>Xr, 1, 0)

}
mean(res)
0.5*(punif(x1, min = 0, max = 2))^2






# Verifico esattezza della forma chiusa di P(x1<Xh<Xr<Y)

B=10^6

res = vector()
x1 = runif(n=1, min=0, max=2)
for(b in 1:B){

  Xh = runif(n=1, min=0, max=2)
  Xr = runif(n=1, min=0, max=2)
  Y = runif(n=1, min=0, max=2)

  res[b] = ifelse(x1<Xh & Xh<Xr & Xr<Y, 1, 0)

}
mean(res)
(1-punif(x1, min = 0, max = 2))^3/6





# Verifico esattezza della forma chiusa di P(Xh<x1<Xr<Y)

B=10^6

res = vector()
x1 = runif(n=1, min=0, max=2)
for(b in 1:B){

  Xh = runif(n=1, min=0, max=2)
  Xr = runif(n=1, min=0, max=2)
  Y = runif(n=1, min=0, max=2)

  res[b] = ifelse(Xh<x1 & x1<Xr & Xr<Y, 1, 0)

}
mean(res)
0.5*(1-punif(x1, min = 0, max = 2))^2*(punif(x1, min = 0, max = 2))





# Verifico esattezza della forma chiusa di P(Xh<Xr<x1<Y)

B=10^6

res = vector()
x1 = runif(n=1, min=0, max=2)
for(b in 1:B){

  Xh = runif(n=1, min=0, max=2)
  Xr = runif(n=1, min=0, max=2)
  Y = runif(n=1, min=0, max=2)

  res[b] = ifelse(Xh<Xr & Xr<x1 & x1<Y, 1, 0)

}
mean(res)
0.5*(1-punif(x1, min = 0, max = 2))*(punif(x1, min = 0, max = 2))^2









# Verifico esattezza della forma chiusa di P(x1<Y, Xh<Y, Xr<Y)

B=10^6

res = vector()
x1 = runif(n=1, min=0, max=2)
for(b in 1:B){

  Xh = runif(n=1, min=0, max=2)
  Xr = runif(n=1, min=0, max=2)
  Y = runif(n=1, min=0, max=2)

  res[b] = ifelse(Xr<Y & Xh<Y & x1<Y, 1, 0)

}
mean(res)
(1-punif(x1, min = 0, max = 2))^3/3+(1-punif(x1, min = 0, max = 2))^2*punif(x1, min = 0, max = 2)+(1-punif(x1, min = 0, max = 2))*(punif(x1, min = 0, max = 2))^2
(1-(punif(x1, min = 0, max = 2))^3)/3




# Verifico esattezza della forma chiusa di
# sum_{h=2}^m (sum_{r=2}^m 1{Xh<Y, Xi<Y, Xr<Y})

B=10^5
m=10

#set.seed(123)
res = vector()

for(b in 1:B){

  Xh = runif(n=(m-1), min=0, max=2)
  Xi = runif(n=1, min=0, max=2)
  X = c(Xi, Xh)
  Y = runif(n=1, min=0, max=2)

  cond1 = ifelse(Xi<Y, 1, 0)
  cond2aus = X<Y
  cond2 = sum(cond2aus%*%t(cond2aus))
  res[b] = cond2*cond1

}

mean(res)
(m-2)*(m/4+1/4)+1/2




# Verifico esattezza della forma chiusa di P(Xr<Xh<Xi<x1<Y)

B=10^6

res = vector()
x1 = runif(n=1, min=0, max=2)
for(b in 1:B){

  Xh = runif(n=1, min=0, max=2)
  Xr = runif(n=1, min=0, max=2)
  Xi = runif(n=1, min=0, max=2)
  Y = runif(n=1, min=0, max=2)

  res[b] = ifelse(Xr<Xh & Xh<Xi & Xi<x1 & x1<Y, 1, 0)

}
mean(res)
(punif(x1, min = 0, max = 2))^3/6*(1-punif(x1, min = 0, max = 2))





# Verifico esattezza della forma chiusa di P(Xr<Xh<x1<Xi<Y)

B=10^6

res = vector()
x1 = runif(n=1, min=0, max=2)
for(b in 1:B){

  Xh = runif(n=1, min=0, max=2)
  Xr = runif(n=1, min=0, max=2)
  Xi = runif(n=1, min=0, max=2)
  Y = runif(n=1, min=0, max=2)

  res[b] = ifelse(Xr<Xh & Xh<x1 & x1<Xi & Xi<Y, 1, 0)

}
mean(res)
(punif(x1, min = 0, max = 2))^2/4*(1-punif(x1, min = 0, max = 2))^2





# Verifico esattezza della forma chiusa di P(Xr<x1<Xh<Xi<Y)

B=10^6

res = vector()
x1 = runif(n=1, min=0, max=2)
for(b in 1:B){

  Xh = runif(n=1, min=0, max=2)
  Xr = runif(n=1, min=0, max=2)
  Xi = runif(n=1, min=0, max=2)
  Y = runif(n=1, min=0, max=2)

  res[b] = ifelse(Xr<x1 & x1<Xh & Xh<Xi & Xi<Y, 1, 0)

}
mean(res)
(punif(x1, min = 0, max = 2))/6*(1-punif(x1, min = 0, max = 2))^3







# Verifico esattezza della forma chiusa di P(x1<Xr<Xh<Xi<Y)

B=10^5

res = vector()
x1 = runif(n=1, min=0, max=2)
for(b in 1:B){

  Xh = runif(n=1, min=0, max=2)
  Xr = runif(n=1, min=0, max=2)
  Xi = runif(n=1, min=0, max=2)
  Y = runif(n=1, min=0, max=2)

  res[b] = ifelse(x1<Xr & Xr<Xh & Xh<Xi & Xi<Y, 1, 0)

}
mean(res)
(1-punif(x1, min = 0, max = 2))^4/24







# Verifico esattezza della forma chiusa di P(x1<Y, Xi<Y, Xh<Y, Xr<Y)

B=10^5

res = vector()
x1 = runif(n=1, min=0, max=2)
for(b in 1:B){

  Xh = runif(n=1, min=0, max=2)
  Xr = runif(n=1, min=0, max=2)
  Xi = runif(n=1, min=0, max=2)
  Y = runif(n=1, min=0, max=2)

  res[b] = ifelse(x1<Y & Xr<Y & Xh<Y & Xi<Y, 1, 0)

}
mean(res)
(1-punif(x1, min = 0, max = 2)^4)/4





# Verifico esattezza della forma chiusa di E[a1^2], i>1 dove
# a1= 1{x1<Y}+sum_{h=1}^m 1{x1<Y, Xh<Y}

B=10^5
m=10

res = vector()
x1 = runif(n=1, min=0, max=2)
for(b in 1:B){

  Xshort = runif(n=(m-1), min=0, max=2)
  X = c(x1, Xshort)
  Y = runif(n=1, min=0, max=2)

  cond2 = X<Y

  # a1
  cond11 = ifelse(x1<Y, 1, 0)
  cond21 = sum(cond2)*cond11
  a1 = cond21+cond11

  res[b] = a1^2
}

mean(res)
4*(1-(punif(x1, min = 0, max = 2)))+
  2*(m-1)*(1-(punif(x1, min = 0, max = 2))^2)+
  (m-1)^2/3*(1-(punif(x1, min = 0, max = 2))^3)





# Verifico esattezza della forma chiusa di E[ai*a1], i>1 dove
# ai= 1{Xi<Y}+sum_{h=1}^m 1{Xi<Y, Xh<Y}
# a1= 1{x1<Y}+sum_{h=1}^m 1{x1<Y, Xh<Y}

B=10^5
m=10

#set.seed(123)
res = vector()
x1 = runif(n=1, min=0, max=2)
for(b in 1:B){

  Xh = runif(n=(m-2), min=0, max=2)
  Xi = runif(n=1, min=0, max=2)
  X = c(x1, Xi, Xh)
  Y = runif(n=1, min=0, max=2)

  cond2 = X<Y

  # a1
  cond11 = ifelse(x1<Y, 1, 0)
  cond21 = sum(cond2)*cond11
  a1 = cond21+cond11

  # ai
  cond1i = ifelse(Xi<Y, 1, 0)
  cond2i = sum(cond2)*cond1i
  ai = cond2i+cond1i

  res[b] = a1*ai

}

mean(res)
(m-2)^2/4*(1-(punif(x1, min = 0, max = 2))^4)+
  2*(m-2)*(1-(punif(x1, min = 0, max = 2))^3)+
  9/2*(1-(punif(x1, min = 0, max = 2))^2)




# Verifico esattezza della forma chiusa di E[ai^2], i>1 dove
# ai= 1{Xi<Y}+sum_{h=1}^m 1{Xi<Y, Xh<Y}

B=10^5
m=10

#set.seed(123)
res = vector()
x1 = runif(n=1, min=0, max=2)
for(b in 1:B){

  Xshort = runif(n=(m-2), min=0, max=2)
  Xi = runif(n=1, min=0, max=2)
  X = c(x1, Xi, Xshort)
  Y = runif(n=1, min=0, max=2)

  cond2 = X<Y

  # ai
  cond1i = ifelse(Xi<Y, 1, 0)
  cond2i = sum(cond2)*cond1i
  ai = cond2i+cond1i

  res[b] = ai^2

}

mean(res)
2+5/3*(m-2)+(m-3)*(m-2)/4+
  2/3*(m-2)*(1-(punif(x1, min = 0, max = 2))^3)+
  5/2*(1-(punif(x1, min = 0, max = 2))^2)



# Verifico esattezza della forma chiusa di
# sum_{h=2}^m (sum_{r=2}^m 1{x1<Y, Xh<Y, Xi<Y, Xr<Y})

B=10^5
m=10

#set.seed(123)
res = vector()
x1 = runif(n=1, min=0, max=2)

for(b in 1:B){

  Xh = runif(n=(m-2), min=0, max=2)
  Xi = runif(n=1, min=0, max=2)
  X = c(Xi, Xh)
  Y = runif(n=1, min=0, max=2)

  cond1 = ifelse(x1<Y & Xi<Y, 1, 0)
  cond2aus = sapply(1:(m-1), function(j) ifelse(X[j]<Y, 1, 0))
  cond2 = sum(cond2aus%*%t(cond2aus))
  res[b] = cond2*cond1

}

mean(res)
(m-2)*(m-3)/4*(1-(punif(x1, min = 0, max = 2))^4)+(m-2)*(1-(punif(x1, min = 0, max = 2))^3)+1/2*(1-(punif(x1, min = 0, max = 2))^2)




# Verifico esattezza della forma chiusa di
# E[(sum_{h=2}^m 1{Xi<Y, Xh<Y})(sum_{r=2}^m 1{Xj<Y, Xr<Y})]

B=10^5
m=10

#set.seed(123)
res = vector()
x1 = runif(n=1, min=0, max=2)

for(b in 1:B){

  Xshort = runif(n=(m-3), min=0, max=2)
  Xi = runif(n=1, min=0, max=2)
  Xj = runif(n=1, min=0, max=2)
  X = c(x1, Xi, Xj, Xshort)
  Y = runif(n=1, min=0, max=2)

  condi = ifelse(Xi<Y, 1, 0)
  condj = ifelse(Xj<Y, 1, 0)
  condX = sum(X<Y)
  condh = condi*condX
  condr = condj*condX
  res[b] = condh*condr

}

mean(res)
1/2*(1-(punif(x1, min = 0, max = 2))^2)+(m-3)/2*(1-(punif(x1, min = 0, max = 2))^4)+4/3*(1-(punif(x1, min = 0, max = 2))^3)+(m-3)*(m-4)/5+5/4*(m-3)+4/3

(12*m^2+21*m+19)/60-(punif(x1, min = 0, max = 2))^2/2-(m-3)/4*(punif(x1, min = 0, max = 2))^4-4/3*(punif(x1, min = 0, max = 2))^3





# Verifico esattezza della forma chiusa di E[ai*aj], i,j>1 i!=j dove
# ai= 1{Xi<Y}+sum_{h=1}^m 1{Xi<Y, Xh<Y}
# aj= 1{Xj<Y}+sum_{r=1}^m 1{Xj<Y, Xr<Y}

B=10^5
m=10

#set.seed(123)
res = vector()
x1 = runif(n=1, min=0, max=2)
for(b in 1:B){

  Xshort = runif(n=(m-3), min=0, max=2)
  Xi = runif(n=1, min=0, max=2)
  Xj = runif(n=1, min=0, max=2)
  X = c(x1, Xi, Xj, Xshort)
  Y = runif(n=1, min=0, max=2)

  cond2 = X<Y

  # aj
  cond1j = ifelse(Xj<Y, 1, 0)
  cond2j = sum(cond2)*cond1j
  aj = cond2j+cond1j

  # ai
  cond1i = ifelse(Xi<Y, 1, 0)
  cond2i = sum(cond2)*cond1i
  ai = cond2i+cond1i

  res[b] = aj*ai

}

mean(res)
3+(m-3)*(4*m+19)/20+1/2*(1-(punif(x1, min = 0, max = 2))^2)+5/3*(1-(punif(x1, min = 0, max = 2))^3)+(m-3)/2*(1-(punif(x1, min = 0, max = 2))^4)





# Verifico esattezza della forma chiusa di E[h(x1, X2,...,Xm,Y)], dove
# h(x1, X2,...,Xm,Y) = sum_{i=1}^m ai dove
# ai= 1{Xi<Y}+sum_{h=1}^m 1{Xi<Y, Xh<Y}

B=10^5
m=10

res = vector()
x1 = runif(n=1, min=0, max=2)
for(b in 1:B){

  Xshort = runif(n=(m-1), min=0, max=2)
  X = c(x1, Xshort)
  Y = runif(n=1, min=0, max=2)

  cond2 = X<Y

  a = vector()
  for(i in 1:m){
    # ai
    cond1i = ifelse(X[i]<Y, 1, 0)
    cond2i = sum(cond2)*cond1i
    a[i] = cond2i+cond1i
  }
  res[b] = sum(a)
}

mean(res)
2*(1-punif(x1, min = 0, max = 2)) +(m-1)*(1-(punif(x1, min = 0, max = 2))^2)+(m-1)/3*(1+m)








# Verifico esattezza della forma chiusa di E[h(x1, X2,...,Xm,Y)^2], dove
# h(x1, X2,...,Xm,Y) = sum_{i=1}^m ai dove
# ai= 1{Xi<Y}+sum_{h=1}^m 1{Xi<Y, Xh<Y}

B=10^5
m=10

res = vector()
x1 = runif(n=1, min=0, max=2)
for(b in 1:B){

  Xshort = runif(n=(m-1), min=0, max=2)
  X = c(x1, Xshort)
  Y = runif(n=1, min=0, max=2)

  cond2 = X<Y

  a = vector()
  for(i in 1:m){
    # ai
    cond1i = ifelse(X[i]<Y, 1, 0)
    cond2i = sum(cond2)*cond1i
    a[i] = cond2i+cond1i
  }
  res[b] = sum(a)^2
}

mean(res)
m^4/5+m^2/2-5/2*m+9/5+4*(1-punif(x1, min = 0, max = 2))+
  (m-1)*(11+m)/2*(1-punif(x1, min = 0, max = 2)^2)+
  (m-1)*(14/3-9)*(1-punif(x1, min = 0, max = 2)^3)+
  ((m-1)*(m-2)*(3/2*m-4)/2)*(1-punif(x1, min = 0, max = 2)^4)







# Verifico esattezza di theta = E[h(X1, X2,...,Xm,Y)]
# dove anche X1 è random

B=10^5
m=10

res = vector()

for(b in 1:B){

  X = runif(n=m, min=0, max=2)
  Y = runif(n=1, min=0, max=2)

  cond2 = X<Y

  a = vector()
  for(i in 1:m){
    # ai
    cond1i = ifelse(X[i]<Y, 1, 0)
    cond2i = sum(cond2)*cond1i
    a[i] = cond2i+cond1i
  }
  res[b] = sum(a)
}

mean(res)
(theta = m+m*(m-1)/3)
mean.incr = sapply(1:B, function(b) mean(res[1:b]))
plot(mean.incr)
abline(h = theta, col = "red")




# Verifico esattezza di E[(h(x1, X2,...,Xm,Y)-theta)^2]

B=10^5
m=10

res = vector()
x1 = runif(n=1, min=0, max=2)
theta = m/3*(m+2)

for(b in 1:B){

  Xshort = runif(n=(m-1), min=0, max=2)
  X = c(x1, Xshort)
  Y = runif(n=1, min=0, max=2)

  cond2 = X<Y

  a = vector()
  for(i in 1:m){
    # ai
    cond1i = ifelse(X[i]<Y, 1, 0)
    cond2i = sum(cond2)*cond1i
    a[i] = cond2i+cond1i
  }
  h = sum(a)
  res[b] = (h - theta)^2
}

mean(res)
4/45*m^4+2/5*m^3+19/30*m^2+47/45*m-7/6+
  (-4/3*m*(m+2)+4)*(1-punif(x1, min = 0, max = 2))+
  (-2/3*m^3-m^2/6+10/3*m-5/2)*(1-(punif(x1, min = 0, max = 2))^2)+
  ((m-1)/3*(8*m-5))*(1-(punif(x1, min = 0, max = 2))^3)+
  (m*(m-2)*(3*m-8))/4*(1-(punif(x1, min = 0, max = 2))^4)
