# ---------------------------- E[aiy^2] ----------------------------

# Verifico esattezza della forma chiusa di E[aiy^2], i>1 dove
# aiy= 1{Xi<y}+sum_{h=1}^m 1{Xi<y, Xh<y}

B=10^5
m=10

y = runif(n=1, min=0, max=2)
res = vector()

for(b in 1:B){

  Xi = runif(n=1, min=0, max=2)
  Xshort = runif(n=(m-1), min=0, max=2)
  X = c(Xi, Xshort)

  cond2 = X<y

  # ai
  cond1i = ifelse(Xi<y, 1, 0)
  cond2i = sum(cond2)*cond1i
  ai = cond2i+cond1i

  res[b] = ai^2

}

mean(res)
4*(punif(y, min = 0, max = 2))+
  5*(m-1)*(punif(y, min = 0, max = 2))^2+
  (m-2)*(m-1)*(punif(y, min = 0, max = 2))^3





# ---------------------------- E[ai_y*aj_y] ----------------------------

# Verifico esattezza della forma chiusa di E[aiy^2], i>1 dove
# aiy= 1{Xi<y}+sum_{h=1}^m 1{Xi<y, Xh<y}

B=10^5
m=10

y = runif(n=1, min=0, max=2)
res = vector()

for(b in 1:B){

  Xi = runif(n=1, min=0, max=2)
  Xj = runif(n=1, min=0, max=2)
  Xshort = runif(n=(m-2), min=0, max=2)
  X = c(Xi, Xj, Xshort)

  cond2 = X<y

  # aj
  cond1j = ifelse(Xj<y, 1, 0)
  cond2j = sum(cond2)*cond1j
  aj = cond2j+cond1j

  # ai
  cond1i = ifelse(Xi<y, 1, 0)
  cond2i = sum(cond2)*cond1i
  ai = cond2i+cond1i

  res[b] = ai*aj

}

mean(res)
9*(punif(y, min = 0, max = 2))^2+
  7*(m-2)*(punif(y, min = 0, max = 2))^3+
  (m-2)*(m-3)*(punif(y, min = 0, max = 2))^4





# ---------------------------- E[h(X1,X2,...,Xm,y)^2] ----------------------------

# Verifico esattezza della forma chiusa di E[h(X1, X2,...,Xm,y)^2], dove
# h(X1, X2,...,Xm,y) = sum_{i=1}^m ai dove
# ai= 1{Xi<y}+sum_{h=1}^m 1{Xi<y, Xh<y}

B=10^5
m=10

res = vector()
y = runif(n=1, min=0, max=2)
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
  res[b] = sum(a)^2
}

mean(res)

4*m*punif(y, min = 0, max = 2)+
  14*(m-1)*m*(punif(y, min = 0, max = 2))^2+
  8*(m-2)*(m-1)*m*(punif(y, min = 0, max = 2))^3+
  m*(m-1)*(m-2)*(m-3)*(punif(y, min = 0, max = 2))^4






# ---------------------------- E[h(X1,X2,...,Xm,y)^2] ----------------------------

# Verifico esattezza della forma chiusa di E[h(X1, X2,...,Xm,y)^2], dove
# h(X1, X2,...,Xm,y) = sum_{i=1}^m ai dove
# ai= 1{Xi<y}+sum_{h=1}^m 1{Xi<y, Xh<y}

B=10^5
m=10

res = vector()
y = runif(n=1, min=0, max=2)
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














