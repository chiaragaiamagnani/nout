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
(closed = 4*(1-(punif(x1, min = 0, max = 2)))+
  5/2*(m-1)*(1-(punif(x1, min = 0, max = 2))^2)+
  (m-1)*(m-2)/3*(1-(punif(x1, min = 0, max = 2))^3))

mean.incr = sapply(1:B, function(b) mean(res[1:b]))
plot(mean.incr)
abline(h = closed, col = "red")


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
(closed = (m-2)^2/4*(1-(punif(x1, min = 0, max = 2))^4)+
  2*(m-2)*(1-(punif(x1, min = 0, max = 2))^3)+
  9/2*(1-(punif(x1, min = 0, max = 2))^2))

(closed2 = (m-2)*(m-3)/4*(1-(punif(x1, min = 0, max = 2))^4)+
    7/3*(m-2)*(1-(punif(x1, min = 0, max = 2))^3)+
    9/2*(1-(punif(x1, min = 0, max = 2))^2))


mean.incr = sapply(1:B, function(b) mean(res[1:b]))
plot(mean.incr)
abline(h = closed, col = "red")
abline(h = closed2, col = "blue")


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
(closed = 2+5/3*(m-2)+(m-3)*(m-2)/4+
  2/3*(m-2)*(1-(punif(x1, min = 0, max = 2))^3)+
  5/2*(1-(punif(x1, min = 0, max = 2))^2))

mean.incr = sapply(1:B, function(b) mean(res[1:b]))
plot(mean.incr)
abline(h = closed, col = "red")



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
(closed = 3+7/4*(m-3)+(m-3)*(m-4)/5+
  7/3*(1-(punif(x1, min = 0, max = 2))^3)+
  (m-3)/2*(1-(punif(x1, min = 0, max = 2))^4))


mean.incr = sapply(1:B, function(b) mean(res[1:b]))
plot(mean.incr)
abline(h = closed, col = "red")


# ---------------------------- theta = E[h(X1,X2,...,Xm,Y)] ----------------------------


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
theta = m+m*(m-1)/3


mean.incr = sapply(1:B, function(b) mean(res[1:b]))
plot(mean.incr)
abline(h = theta, col = "red")

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
(closed = 2*(1-punif(x1, min = 0, max = 2))+
  (m-1)*(1-(punif(x1, min = 0, max = 2))^2)+
  (m-1)/3*(1+m))

mean.incr = sapply(1:B, function(b) mean(res[1:b]))
plot(mean.incr)
abline(h = closed, col = "red")



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

# closed sbagliato, giusto closed2
(closed = (m-1)*(m-2)/15*(3*m^2+9*m+16)+
  4*(1-punif(x1, min = 0, max = 2))+
  14*(m-1)*(1-punif(x1, min = 0, max = 2)^2)+
  8*(m-1)*(m-2)*(1-punif(x1, min = 0, max = 2)^3)+
  (m-1)*(m-2)*(m-3)*(1-punif(x1, min = 0, max = 2)^4))

(closed2 = 5/6-4*punif(x1, min = 0, max = 2)+
  5*punif(x1, min = 0, max = 2)^2-2*punif(x1, min = 0, max = 2)^3+
  7*m/4-5*punif(x1, min = 0, max = 2)^2*m+
  3*punif(x1, min = 0, max = 2)^3*m+7/6*m^2-punif(x1, min = 0, max = 2)^3*m^2+m^3/4+
  1/20*(m-1)*(14-180*punif(x1, min = 0, max = 2)^2-140*punif(x1, min = 0, max = 2)^3*(m-2)+29*m+
                19*m^2+4*m^3-20*punif(x1, min = 0, max = 2)^4*(6-5*m+m^2)))

(closed3 = (m-1)/15*(3*m^3+3*m^2-2*m-2)+
    4*(1-punif(x1, min = 0, max = 2))+
    14*(m-1)*(1-punif(x1, min = 0, max = 2)^2)+
    8*(m-1)*(m-2)*(1-punif(x1, min = 0, max = 2)^3)+
    (m-1)*(m-2)*(m-3)*(1-punif(x1, min = 0, max = 2)^4))


mean.incr = sapply(1:B, function(b) mean(res[1:b]))
plot(mean.incr)
abline(h = closed, col = "red")
abline(h = closed2, col = "blue")
abline(h = closed3, col = "darkgreen")


# -------------------------- z11 ----------------------------------------------

B=10^5
m = 10

theta = m+m*(m-1)/3
theta2 = theta^2

Eh2x1 = function(x1){
  out = (m-1)/15*(3*m^3+3*m^2-2*m-2)+
    4*(1-punif(x1, min = 0, max = 2))+
    14*(m-1)*(1-punif(x1, min = 0, max = 2)^2)+
    8*(m-1)*(m-2)*(1-punif(x1, min = 0, max = 2)^3)+
    (m-1)*(m-2)*(m-3)*(1-punif(x1, min = 0, max = 2)^4)
  return(out)
}

Ehx1 = function(x1){
  out = 2*(1-punif(x1, min = 0, max = 2))+
    (m-1)*(1-(punif(x1, min = 0, max = 2))^2)+
    (m-1)/3*(1+m)
  return(out)
}

x1s = runif(n=B, min=0, max=2)
res = sapply(x1s, function(x) Eh2x1(x)+theta2-2*theta*Ehx1(x))


mean(res)

(z11 = 4/45*m^4+16/45*m^3+19/45*m^2+2/15*m)

mean.incr = sapply(1:B, function(b) mean(res[1:b]))
plot(mean.incr)
abline(h = z11, col = "red")



























