
# ------------------------------------------------------------------------ #
#  Implementation of generalized WMW tests with asymptotical distribution  #
# ------------------------------------------------------------------------ #



#' stat.Tk
#'
#' @description It computes the ranks of the test observations in the pooled
#' vector of calibration and test scores. Then, component-wisely compute the power
#' from \eqn{1} to \eqn{k} of the rank vector and component-wisely sum them
#'
#' @param Z : pooled score vector with the first \eqn{m} components corresponding
#' to calibration observations and the last \eqn{n} components corresponding
#' to test observations
#' @param m : calibration sample size
#' @param k : order of the LMPI test statistic
#'
#' @return Given the pooled score vector \eqn{Z=(X,Y)} where \eqn{X} is the calibration score
#' vector and \eqn{Y} is the test score vector, for each observation in the test sample
#' it returns
#' \deqn{T_{k,i}=\sum{j=1}^m\mathbb{1}\{X_j<Y_i\}+\left(\sum{j=1}^m\mathbb{1}\{X_j<Y_i\}\right)^2+\ldots+\left(\sum{j=1}^m\mathbb{1}\{X_j<Y_i\}\right)^k.}
#'
#'
stat.Tk <- function(Z, m, k) {
  m = as.double(m)
  k = as.double(k)
  N = as.double(length(Z))
  n = as.double(N-m)
  R = rank(Z)[(m+1):N]-1

  if(k>1){ # test statistic Shiraishi (1985)
    if(n==1){
      Tk_i = (k+1)/((N+1):(N+k))*prod(sapply(1:k, function(l){R+l-1}))
    } else {
      Tk_i = (k+1)/((N+1):(N+k))*apply(sapply(1:k, function(l){R+l-1}), MARGIN=1, FUN = prod)
    }
  } else if(k==1){ # classic rank sum Wilcoxon test statistic
    Tk_i = R
  }
  return(Tk_i)
}






# Compute the k-th moment of a Beta(h, N-h+1), k =1,2,3,...
#' beta_moment.k
#'
#' @param h : an integer between \eqn{1} and \eqn{N}
#' @param N : integer
#' @param k : order of moment
#'
#' @return \eqn{k}the moment of a Beta(h, N-h+1) distribution
#' @export
#'
#' @examples
#' mom = beta_moment.k(h=5,N=10,k=3)
beta_moment.k = function(h,N,k){

  stopifnot(k>=1 & k%%1==0)

  r.seq = seq(from=0,to=k-1)
  num = h+r.seq
  den = N+1+r.seq
  out = prod(num/den)

  return(out)
}



#' mean_analytical.Lk
#'
#' @param N : pooled sample size
#' @param k : an integer greater than 1, order of the Lehmann's alternative distribution
#' @param n : test sample size
#'
#' @return Asymptotic mean of Shiraishi (1985) test statistic
#' @export
#'
#' @examples
#' m = mean_analytical.Lk(N=100,k=2,n=50)
mean_analytical.Lk = function(N,k,n){

  stopifnot(k>=2 & k%%1==0)

  aN_j = k*mapply(1:N, FUN=function(h) beta_moment.k(h,N=N,k=k-1))
  out = n*mean(aN_j)

  return(out)

}



#' var_analytical.Lk
#'
#' @param N : pooled sample size
#' @param k : an integer greater than 1, order of the Lehmann's alternative distribution
#' @param n : test sample size
#'
#' @return Asymptotic variance of Shiraishi (1985) test statistic
#' @export
#'
#' @examples
#' v = var_analytical.Lk(N=100,k=2,n=50)
var_analytical.Lk = function(N,k,n){

  stopifnot(k>=2 & k%%1==0)

  m=N-n

  aN_j = k*mapply(1:N, FUN=function(h) beta_moment.k(h,N=N,k=k-1))
  bar.aN = mean(aN_j)
  out = n*m*sum((aN_j-bar.aN)^2)/(N*(N-1))

  return(out)

}







#' asymptotic.moments.Tk
#'
#' @description It computes the mean and the variance of the asymptotic distribution of \eqn{T_k}.
#' For \eqn{k=1,2,3} the exact mean is given, i.e., the finite-sample mean,
#' and the provided estimate of the variance is more accurate than the first order asymptotic approximation.
#'
#' @param m : calibration sample size
#' @param n : test sample size
#' @param k : order of the LMPI test statistic
#'
#' @return It returns the mean and the variance of the asymptotic distribution of \eqn{T_k}.
#' For \eqn{k=1,2,3} the exact mean is given, i.e., the finite-sample mean,
#' and the provided estimate of the variance is more accurate than the first order asymptotic approximation.
asymptotic.moments.Tk <- function(m, n, k) {

  stopifnot(k>=1 & k%%1==0)

  N = as.double(n+m)
  if(k==1){

    mean.Tk = m*n/2 + n*(n-1)/2
    variance.Tk = m*n*(m+n)/12

  } else if(k>1){
    # Asymptotic distribution Shiraishi (1985)
    mean.Tk = mean_analytical.Lk(N=N,k=k,n=n)
    variance.Tk = var_analytical.Lk(N=N,k=k,n=n)
  }

  moments.Tk = list("mean.Tk"=mean.Tk, "variance.Tk"=variance.Tk)

  return(moments.Tk)
}





#' asymptotic.critical.Tk
#'
#' @description It computes the \eqn{(1-\alpha)}-quantile of LMPI \eqn{T_k} test statistic
#' based on asymptotic normal approximation. For \eqn{k=1,2,3} the exact mean and a more accurate estimate are used.

#' @param m : calibration sample size
#' @param n : test sample size
#' @param k : order of the LMPI test statistic
#' @param alpha : significance level. Default value is set equal to 0.1
#'
#'
#' @return It returns the \eqn{(1-\alpha)}-quantile of LMPI \eqn{T_k} test statistic
#' based on asymptotic normal approximation.
#'
asymptotic.critical.Tk <- function(m, n, k, alpha=0.1) {

  stopifnot(k>=1 & k%%1==0)

  moments.Tk = asymptotic.moments.Tk(m=m,n=n,k=k)
  mean.Tk = as.double(moments.Tk$mean.Tk)
  variance.Tk = as.double(moments.Tk$variance.Tk)

  critical.value = as.double(stats::qnorm(alpha, mean=mean.Tk, sd = sqrt(variance.Tk), lower.tail = F))

  return(critical.value)
}



#' asymptotic.pvalue.Tk
#'
#' @description It computes the approximated *p*-value of the LMPI \eqn{T_k}
#' test statistic based on the asymptotic normal approximation.
#'
#' @param m : calibration sample size
#' @param n : test sample size
#' @param k : order of the LMPI test statistic
#' @param T.obs : observed value of the test statistic
#'
#'
#' @return It returns the approximated *p*-value of the LMPI \eqn{T_k}
#' test statistic based on the asymptotic normal approximation
#'
asymptotic.pvalue.Tk <- function(m, n, k, T.obs) {

  stopifnot(k>=1 & k%%1==0)

  moments.Tk = asymptotic.moments.Tk(m=m,n=n,k=k)
  mean.Tk = moments.Tk$mean.Tk
  variance.Tk = moments.Tk$variance.Tk

  p.value = stats::pnorm(q=T.obs, mean=mean.Tk, sd = sqrt(variance.Tk), lower.tail = F)

  return(p.value)
}


