
#' d_ttest
#'
#' @param S_Y : score vector of test observations
#' @param S_X : score vector of calibration observations
#' @param alpha : significance level of the local test. Default value is set equal to 0.1.
#'
#' @return An integer which is the \eqn{(1 − \alpha)}-confidence lower bound for the
#' number of true discoveries using step-down Dunnett test procedure, which is
#' closed testing procedure using Dunnett local test.
#' The selection set, i.e. the set of hypothesis indices that we are
#' interested in is \eqn{[m]=:\{1,...,m\}} by default.
#' @export
#'
#' @examples
#' #' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' d_Simes(S_Y=Sy, S_X=Sx)
#'
#'
d_ttest = function(S_Y, S_X, alpha = 0.1){
  m = length(S_Y)
  n = length(S_X)
  z = c(S_X,S_Y)

  group = as.factor(c(rep(0,m),1:n))
  dat = data.frame(z,group)
  fit = stats::aov(z ~ group, data = dat)
  dun_onesided = multcomp::glht(fit,
                                linfct = multcomp::mcp(group = "Dunnett"),
                                alternative = "greater")
  sum.free = summary(dun_onesided, test = multcomp::adjusted(type = "free"))
  d = sum(sum.free$test$pvalues<=alpha)

  return(d)
}

