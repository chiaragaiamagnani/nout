#' crit.WMW
#'
#' @description Given \eqn{m} observations in the calibration set and \eqn{n} observations in the test set,
#' it relies on R function ***qwilcox*** and returns the vector of critical values of the
#' Wilcoxon-Mann-Whitney test statistic for each level of the closed testing.
#'
#'
#' @param m : number of observations in the calibration set
#' @param n : number of observations in the test set
#' @param alpha : significance level of the local test. Default value is set equal to 0.1
#' @param n.exact : maximum value of the sample size of the second sample for which the critical
#' values of the Wilcoxon-Mann-Whitney statistic are exactly computed using ***qwilcox*** function.
#' Default value is set equal to 10.
#' @param exact : logical value. If TRUE, exact computation of critical values is performed for
#' all values of \eqn{k\in{1,\ldots,m}} when \eqn{n\leq 20}. Otherwise, exact computation of
#' critical values is performed for all values of \eqn{k\in\{1,\ldots,n.exact\}} when \eqn{n<1000}.
#' Default value is FALSE.
#'
#'
#'
#' @return A R object of class *crit.val.info*, which is a list consisting
#' of the following elements: \itemize{
#' \item m number of test observations
#' \item n number of calibration observations
#' \item alpha significance level of the local test
#' \item crit.val vector of critical values of the Wilcoxon-Mann-Whitney
#' statistic \eqn{U(n,k)}, \eqn{k=m,\ldots,1} (first element of the vector
#' corresponds to \eqn{k=m} and the last one to \eqn{k=1.} Given two samples
#' \eqn{(X_1,\ldots,X_n)} and \eqn{(Y_1,\ldots,Y_k)}, the Wilcoxon-Mann-Whitney
#' statistic is computed as \deqn{U(n,k) = \sum_{i = 1}^{n}\sum_{j = 1}^{k} \mathbb{1}\{X_i > Y_j\}}
#' and it has mean \eqn{kn/2} and variance \eqn{kn(k+n+1)/12}.
#'
#' When \eqn{n<1000} and \eqn{m<10} the value of the Wilcoxon-Mann-Whitney statistic
#' is exactly computed with the R function ***qwilcox***; otherwise, we use
#' the normal approximation using continuity correction \deqn{U(n,k) \sim N(kn/2-1/2, kn(k+n+1)/12)}.
#' }
#'
#'
#' @export
#'
#' @examples
#' crit.WMW(n=100,m=2000,alpha=0.1)
#'
#' crit.WMW(n=100,m=900,alpha=0.1)
#'
#' n=10;m=20
#' crits = crit.WMW(m=m, n=n,alpha=0.1)$crit.val
#' plot(x=1:n, y=crits, xlab="n-k+1", main = "Critical values of WMW statistic for m=20 and k=1,...,n")
#'
#'
#'
crit.WMW = function(m, n, alpha, n.exact=10, exact=F){

  # exact=F
  if(exact==F){
    if(m>=10^3){
      crit = sapply(1:n, function(k) stats::qnorm(alpha, mean=((n-k+1)*m/2-0.5),
                                                  sd = sqrt((n-k+1)*m*((n-k+1)+m+1)/12),
                                                  lower.tail = F))
    }
    else{
      if(n>n.exact){
        nn=n-n.exact
        crit2 = sapply(1:nn, function(k) stats::qnorm(alpha, mean=((n-k+1)*m/2-0.5),
                                                      sd = sqrt((n-k+1)*m*((n-k+1)+m+1)/12),
                                                      lower.tail = F))

        crit1 = sapply((nn+1):n, function(k) stats::qwilcox(p=alpha, n=n-k+1, m=m,
                                                            lower.tail = FALSE))
        crit=c(crit2, crit1)
      }
      else{
        crit = sapply(1:n, function(k) stats::qwilcox(p=alpha, n=n-k+1, m=m,
                                                      lower.tail = FALSE))
      }
    }
    res = list("m" = m, "n" = n, "crit.val" = crit, "alpha" = alpha)
    class(res) = "crit.val.info"
  }

  # exact=T
  else{
    if(m>=10^3){
      crit = sapply(1:n, function(k) stats::qnorm(alpha, mean=((n-k+1)*m/2-0.5),
                                                  sd = sqrt((n-k+1)*m*((n-k+1)+m+1)/12),
                                                  lower.tail = F))
    }

    if(m>20 & m<10^3){
      if(n>n.exact){
        nn=n-n.exact
        crit2 = sapply(1:nn, function(k) stats::qnorm(alpha, mean=((n-k+1)*m/2-0.5),
                                                      sd = sqrt((n-k+1)*m*((n-k+1)+m+1)/12),
                                                      lower.tail = F))

        crit1 = sapply((nn+1):n, function(k) stats::qwilcox(p=alpha, n=n-k+1, m=m,
                                                            lower.tail = FALSE))
        crit=c(crit2, crit1)
      }
      else{
        crit = sapply(1:n, function(k) stats::qwilcox(p=alpha, n=n-k+1, m=m,
                                                      lower.tail = FALSE))
      }
    }

    if(m<=20){
      crit = sapply(1:n, function(k) stats::qwilcox(p=alpha, n=n-k+1, m=m,
                                                    lower.tail = FALSE))
    }
  }

  res = list("m" = m, "n" = n, "crit.val" = crit, "alpha" = alpha)
  class(res) = "crit.val.info"

  return(res)
}






#' exact.crit.WMW
#'
#' @description Given \eqn{m} observations in the calibration set and \eqn{n} observations in the test set,
#' it relies on R function ***qwilcox*** and returns the vector of critical values of the
#' Wilcoxon-Mann-Whitney test statistic for each level of the closed testing.
#'
#'
#' @param m : number of observations in the calibration set
#' @param n : number of observations in the test set
#' @param alpha : significance level of the local test. Default value is set equal to 0.1
#'
#'
#' @return A R object of class *crit.val.info*, which is a list consisting
#' of the following elements: \itemize{
#' \item m number of calibration observations
#' \item n number of test observations
#' \item alpha significance level of the local test
#' \item crit.val vector of critical values of the Wilcoxon-Mann-Whitney
#' statistic \eqn{U(n,k)}, \eqn{k=m,\ldots,1} (first element of the vector
#' corresponds to \eqn{k=m} and the last one to \eqn{k=1.} Given two samples
#' \eqn{(X_1,\ldots,X_n)} and \eqn{(Y_1,\ldots,Y_k)}, the Wilcoxon-Mann-Whitney
#' statistic is computed as \deqn{U(n,k) = \sum_{i = 1}^{n}\sum_{j = 1}^{k} \mathbb{1}\{X_i > Y_j\}}
#' and it has mean \eqn{kn/2} and variance \eqn{kn(k+n+1)/12}.}
#'
#'
#'
#' @export
#'
#' @examples
#'
#' exact.crit.WMW(n=10,m=5,alpha=0.1)
#'
#' n=8;m=10
#' crits = exact.crit.WMW(m=m, n=n,alpha=0.1)$crit.val
#' plot(x=1:n, y=crits, xlab="n-k+1", main = "Critical values of WMW for m=10 and k=1,...,n=8")
#'
#'
#'
exact.crit.WMW = function(m, n, alpha=0.1){

  crit = sapply(1:n, function(k) stats::qwilcox(p=alpha, n=n-k+1, m=m,
                                                lower.tail = FALSE))
  res = list("m" = m, "n" = n, "crit.val" = crit, "alpha" = alpha)
  class(res) = "crit.val.info"

  return(res)
}







#' d_t2
#'
#' @description It returns the lower bound for the number of true discoveries in closed testing procedure
#' using Wilcoxon-Mann-Whitney local test.
#'
#' @param S_Y : score vector for the test set
#' @param S_X : score vector for the calibration set
#' @param n.exact : maximum value of the sample size of the second sample for which the critical
#' values of the Wilcoxon-Mann-Whitney statistic are exactly computed using ***qwilcox*** function.
#' Default value is set equal to 10.
#' @param alpha : significance level. Default level is set equal to 0.1
#'
#' @return An integer which is the \eqn{(1 − \alpha)}-confidence lower bound for
#' the number of true discoveries in closed testing procedure using
#' Wilcoxon-Mann-Whitney local test that rejects the global null hypothesis at level \eqn{\alpha}
#' when the test statistic \deqn{\sum_{i=1}^n\sum_{j=1}^m \mathbb{1}\{X_j<Y_i\}} is greater than
#' the critical value corresponding to the significance level \eqn{\alpha}.
#' The selection set is trivial, i.e., we are interested in testing all the observations in the test set by default.
#'
#'
#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70)
#' Sy = setdiff(Sxy, Sx)
#' crit = crit.WMW(m=length(Sx), n=length(Sy), alpha=0.1)
#' d_t2(S_Y=Sy, S_X=Sx, alpha=0.1)
#'
d_t2 = function(S_Y,S_X,n.exact=10,alpha=0.1){

  # if (!inherits(crit, "crit.val.info")) {
  #   stop("Error: crit class not correct")
  # }

  # alpha = crit$alpha
  # critical_vals = crit$crit.val
  # m = crit$m
  # n = crit$n

  # if(length(S_Y)!=n){
  #   stop("Error: length of S_Y differs from m")
  # }
  #
  # if(length(S_X)!=m){
  #   stop("Error: length of S_X differs from n")
  # }

  n = length(S_Y)
  m = length(S_X)

  if(min(n,m)<=n.exact){
    if(n.exact>20){
      stop("Execution stopped: exact computation of critical values is computationally too demanding.")
    }
    crit = nout::crit.WMW(m=m, n=n, n.exact=n.exact, alpha=alpha)$crit.val
  }

  if(min(n,m)>n.exact){
    crit = nout::crit.WMW(m=m, n=n, n.exact=n.exact, alpha=alpha)$crit.val
  }
  # Ranks of S_Y[i] in (S_X,S_Y[i])
  U_i = sort(sapply(1:n, function(i) sum(S_Y[i]>S_X)),
             decreasing = TRUE)

  # For each k in {1,...,m} consider the worst case scenario
  # (consider the k smallest U_i)
  U = sapply(1:n, function(k) sum(U_i[k:n]))

  d = sum(cumsum(U >= crit) == 1:n)

  return(d)
}



