# Density function of some of the Lehmann's alternatives
## WMW k=2 (classic WMW)
g1 = function(x, k=2){
  out = ifelse(x<1 & x>0, k*x^(k-1), 0)
  return(out)
}

rg1 = function(n, rg_null, k=2){
  out = replicate(n, max(rg_null(k)))
  return(out)
}

## WMW k=3
g2 = function(x, k=3){
  out = ifelse(x<1 & x>0, k*x^(k-1), 0)
  return(out)
}

rg2 = function(n, rg_null, k=3){
  out = replicate(n, max(rg_null(k)))
  return(out)
}


## WMW k=10
g9 = function(x, k=10){
  out = ifelse(x<1 & x>0, k*x^(k-1), 0)
  return(out)
}

rg9 = function(n, rg_null, k=10){
  out = replicate(n, max(rg_null(k)))
  return(out)
}

