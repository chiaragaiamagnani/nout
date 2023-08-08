
findhStoreySimes = function(pvalues, lambda = 0.5, alpha){

  n = length(pvalues)
  p.ordered = sort(pvalues, decreasing=F)
  # pi0Storey.mod  vector of length n where the first one corresponds to s=n and the last one to s=1
  # pi0Storey.mod = sapply(n:1, function(s) (1+sum(p.ordered[(n-s+1):n]>=lambda))/(s*(1-lambda)))

  # Storey estimator will be used in the closed testing procedure
  # in every levels except for lowest ones, when the set of
  # considered pvalues has cardinality less than or equal to 2.
  pi0Storey.mod.highlevels = sapply(n:3, function(s) (1+sum(p.ordered[(n-s+1):n]>=lambda))/(s*(1-lambda)))
  pi0Storey.mod = c(pi0Storey.mod.highlevels, 1,1)
  # coeffs =  vector of length n where the first one corresponds to s=n and the last one to s=1
  coeffs = alpha/((n:1)*pi0Storey.mod)

  h=0; cont=T; s=n;

  while(cont==T & s>0){
    p = p.ordered[((n-s+1)):n]
    # thr = vector of length s where the first one has weight 1 and the last one has weight equal to s
    thr = (1:s)*coeffs[(n-s+1)]
    if(sum(p>thr)==s){
      h=s; cont=F
    }
    else
      s=s-1
  }

  if(h==0){
    h.res = list("h" = h,
                 "n" = n,
                 "pvalues" = pvalues,
                 "coeff" = alpha*(1-lambda))
  }
  else{
    h.res = list("h" = h,
                 "n" = n,
                 "pvalues" = pvalues,
                 # "coeff" = coeffs[which(coeffs == h)]
                 "coeff" = coeffs[n-h+1])
  }

  class(h.res) = "hres.class"

  return(h.res)

}



finddStoreySimes = function(input, S, alpha){

  if (!inherits(input, "hres.class")) {
    stop("Error: input class must be hres.class")
  }

  n = input$n
  p.selected = input$pvalues[S]
  coeff = input$coeff

  s = length(S)
  ds = vector()

  for(u in 1:s){
    thr = u*coeff
    ds[u] = 1-u+sum(p.selected<=thr)
  }

  return(max(ds))
}






#' dselection_Simes
#'
#' @param S_X vector of calibration scores
#' @param S_Y vector of test scores
#' @param S set of indices corresponding to the *p*-values to be tested
#' @param lambda : parameter involved in the computation of Storey estimator.
#' Default value is set equal to 0.5
#' @param alpha significance level
#'
#' @return An integer which is the \eqn{(1 − \alpha)}-confidence lower bound for
#' the number of true discoveries in closed testing procedure using
#' Simes local test with Storey's estimator for the proportion of true null
#' hypotheses applied to conformal *p*-values. The selection set, i.e. the set of hypothesis
#' indices that we are interested in can be chosen.
#'
#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' dselection_StoreySimes(S_Y=Sy, S_X=Sx, S = 1:29, alpha = 0.6)
#'
dselection_StoreySimes = function(S_X, S_Y, S, lambda = 0.5, alpha){

  n = length(S_Y)
  m = length(S_X)
  pval = sapply(1:n, function(i) (1+sum(S_X >= S_Y[i]))/(m+1))

  h.res = findhStoreySimes(pvalues = pval, alpha = alpha)
  d = finddStoreySimes(input = h.res, S=S, alpha = alpha)

  return(d)
}

