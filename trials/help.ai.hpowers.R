# ---------------------------- E[a1^2] ----------------------------

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
  5/2*(m-1)*(1-(punif(x1, min = 0, max = 2))^2)+
  (m-1)*(m-2)/3*(1-(punif(x1, min = 0, max = 2))^3)





# ---------------------------- E[a1*ai] ----------------------------

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





# ---------------------------- E[ai^2] ----------------------------

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



# ---------------------------- E[ai*aj] ----------------------------

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
3+7/4*(m-3)+(m-3)*(m-4)/5+
  7/3*(1-(punif(x1, min = 0, max = 2))^2)+
  (m-3)/2*(1-(punif(x1, min = 0, max = 2))^4)




# ---------------------------- E[h(X1,X2,...,Xm,Y)] ----------------------------


# Verifico esattezza della forma chiusa di E[h(X1, X2,...,Xm,Y)], dove
# anche X1 è random e
# h(x1, X2,...,Xm,Y) = sum_{i=1}^m ai dove
# ai= 1{Xi<Y}+sum_{h=1}^m 1{Xi<Y, Xh<Y}

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
m+m*(m-1)/3



# ---------------------------- E[h(x1,X2,...,Xm,Y)] ----------------------------

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
2*(1-punif(x1, min = 0, max = 2))+
  (m-1)*(1-(punif(x1, min = 0, max = 2))^2)+
  (m-1)/3*(1+m)





# ---------------------------- E[h(x1,X2,...,Xm,Y)^2] ----------------------------

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

1/15*(m^2-1)*(3*m^2-2)+
  4*(1-punif(x1, min = 0, max = 2))+
  (m-1)*14*(1-punif(x1, min = 0, max = 2)^2)+
  (m-1)*(m-2)*22/3*(1-punif(x1, min = 0, max = 2)^3)+
  (m-1)*(m-2)*(2*m-5)/2*(1-punif(x1, min = 0, max = 2)^4)




# ---------------------------- E[h(X1,X2,...,Xm,y)] ----------------------------

# Verifico esattezza della forma chiusa di E[h(X1, X2,...,Xm,y)], dove
# h(X1, X2,...,Xm,y) = sum_{i=1}^m ai dove
# ai= 1{Xi<y}+sum_{h=1}^m 1{Xi<y, Xh<y}

B=10^5
m=10

y = runif(n=1, min=0, max=2)

res = vector()
for(b in 1:B){

  X = runif(n=m, min=0, max=2)

  cond2 = X<y

  a = vector()
  for(i in 1:m){
    # ai
    cond1i = ifelse(X[i]<y, 1, 0)
    cond2i = sum(cond2)*cond1i
    a[i] = cond2i+cond1i
  }
  res[b] = sum(a)
}

mean(res)

2*m*punif(y, min = 0, max = 2)+
  (m-1)*m*(punif(y, min = 0, max = 2))^2







