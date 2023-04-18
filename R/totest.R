
totest = function(S_X, S_Y, lambda=0.5, alpha=0.1){

  m = length(S_Y)
  n = length(S_X)
  pval = sort(sapply(1:m, function(i) (1+sum(S_X >= S_Y[i]))/(n+1)), decreasing = F)

  # test d_Simes
  test.Simes=T
  k=0
  while(test.Simes==T & k<m){
    k=k+1
    coeff = seq(from=1, to=m-k+1, by=1)
    pval.Simes = min(pval[k:m]/coeff)
    thr = alpha/(m-k+1)
    test.Simes = pval.Simes<=thr
  }
  d_Simes = k-1


  # test d_StoreySimes
  test.StoSimes=T
  k=0
  while(test.StoSimes==T & k<m){
    k=k+1
    coeff = seq(from=1, to=k, by=1)
    pval.Simes = min(pval[k:m]/coeff)
    pi0 = (1+sum(pval[k:m]>lambda))/((m-k+1)*(1-lambda))
    thr = alpha/((m-k+1)*pi0)
    test.StoSimes = pval.Simes<=thr
  }
  d_StoSimes = k-1

  return(c("d_Simes" = d_Simes, "d_StoSimes"=d_StoSimes))
}


S_X=c(10,9,8,2,12,7,14,5,13,4)
S_Y = c(1,3,6)
(true.val = totest(S_X, S_Y))
d_Simes(S_Y, S_X)
d_StoreySimes(S_Y, S_X)
crit = compute_critWMW(m=length(S_Y),n=length(S_X),alpha=0.1)
d_mannwhitney(S_Y,S_X,crit)



S_X=c(10,9,8,17,12,7,14,5,13,4)
S_Y = c(1,3,6,2)
(true.val = totest(S_X, S_Y))
d_Simes(S_Y, S_X)
d_StoreySimes(S_Y, S_X)
crit = compute_critWMW(m=length(S_Y),n=length(S_X),alpha=0.1)
d_mannwhitney (S_Y,S_X,crit)





S_X=c(10,9,8,17,20,12,7,14,5,13,4,16,18,24,100,25,23,22,30,31,33,28,29)
S_Y = c(1,3,6,2)
(true.val = totest(S_X, S_Y))
d_Simes(S_Y, S_X)
d_StoreySimes(S_Y, S_X)
crit = compute_critWMW(m=length(S_Y),n=length(S_X),alpha=0.1)
d_mannwhitney(S_Y,S_X,crit)


