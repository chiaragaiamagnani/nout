test_that("find_d works", {
  g2 = function(x, k=2) ifelse(x<1 & x>0, k*x^(k-1), 0)
  rg2 = function(rnull, k=2) max(rnull(k))
  X = runif(50)
  Y = replicate(50, rg2(rnull=runif))
  res = find_d(X, Y, B=100)
  res = find_d(X, Y, local.test="higher", k=3, B=100)
  res = find_d(X, Y, local.test="g", g.hat = g2, S=c(1:20), monotonicity=TRUE, B=100)

})
