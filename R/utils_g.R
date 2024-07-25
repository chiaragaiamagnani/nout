
# -------------------------------------------------------------------------------------------- #
#  Estimate the mixing proportion in the mixture model using the method in Petra & Sen (2016)  #
# -------------------------------------------------------------------------------------------- #

# This function calculates the distance $\gamma  \ d_n(\hat{F}_{s,n}^{\gamma},\check{F}_{s,n}^\gamma)$
# for grid of gamma values in [0,1].
# data is a numeric vector containing observations from the mixture model.
# Bigger gridsize  gives more  accurate estimates of alpha.

EstMixMdl <- function(data,Fb,gridsize=200)
{
  n <- length(data)    		 			 ## Length of the data set
  data <- sort(data)       ## Sorts the data set
  data.1 <- unique(data)	 	## Finds the unique data points
  Fn <- stats::ecdf(data)			        ## Computes the empirical DF of the data
  Fn.1 <- Fn(data.1)		      ## Empirical DF of the data at the data points
  Fb.data <- Fb(data.1)       ## evaluating Fb at the the unique data points
  ## Compute the weights (= frequency/n) of the unique data values, i.e., dF_n
  Freq <- diff(c(0,Fn.1))
  distance <- rep(0,gridsize)
  distance[0]<- sqrt(t((Fn.1-Fb.data)^2)%*%Freq)
  for(i in 1:gridsize)
  {
    a <- i/gridsize               ## Assumes a value of the mixing proportion
    F.hat <- (Fn.1-(1-a)*Fb.data)/a			    ## Computes the naive estimator of F_s
    F.is <- Iso::pava(F.hat,Freq,decreasing=FALSE)	## Computes the Isotonic Estimator of F_s
    F.is[which(F.is<=0)] <- 0
    F.is[which(F.is>=1)] <- 1
    distance[i] <- a*sqrt(t((F.hat-F.is)^2)%*%Freq);
  }
  return(distance)
}



# The following function evaluates the numerical second derivative of any function
Comp_2ndDer <- function(dist.alpha, gridsize)
{
  dder <- diff(dist.alpha)    ## Computes the 1st order differences
  dder <- diff(dder)      ## Computes the 2nd order differences
  dder <- c(0,0,dder)       ## The numerical double derivative vector

  return(dder)
}


estimate_mixing_prop = function(X, Y, F_null, gridsize=4000){

  n <- length(Y)
  # m <- length(X)
  # pval = sapply(1:n, function(i) (1+sum(X >= Y[i]))/(m+1))
  # dist.alpha <- EstMixMdl(Y,F_null,gridsize)
  c.n<-0.1*log(log(n))
  # Est<- sum(dist.alpha>c.n/sqrt(n))/gridsize

  dist.out = mixmodel::dist.calc(data = Y, gridsize = 2000)
  alp.hat <- sum(dist.out$distance > c.n/ sqrt(n))/gridsize
  return(alp.hat)

}



# -------------------------------------------------------------------------------- #
#                                    Estimate g                                    #
# -------------------------------------------------------------------------------- #
compute_estimate_null_distr = function(X){

  F.hat = stats::ecdf(X)
  f.hat = statip::densityfun(X, kernel="rectangular")

  return(list("ecdf" = F.hat, "edensity"=f.hat))
}


KDE_mixture_density = function(X,Y,null_cdf,ker){

  stopifnot("Error: kernel must be in .kernelsList()"= ker%in%statip::.kernelsList())

  # Estimate mixture density
  Z = c(X,Y)
  FZ = null_cdf(Z)
  mixture_hat = statip::densityfun(FZ, kernel=ker)

  return(mixture_hat)
}




invert_mixture = function(mixture_density, null_density, prop.out){

  g.hat = function(x){

    # Invert and find estimated outlier distribution density
    stopifnot("Error: x must be in [0,1]"= 0<=x & x<=1)

    g.hat_eval = (mixture_density(x)-prop.out*null_density(x))/(1-prop.out)

    res = ifelse(g.hat_eval>0,g.hat_eval,0)
    return(res)
  }

  return(g.hat)
}



estimate_g = function(X1,X2,Y, constraint=NULL, ker="uniform"){

  if(is.null(constraint)){
    stopifnot("Error: kernel must be in .kernelsList()"= ker%in%statip::.kernelsList())

    # F.hat = compute_estimate_null_distr(X=X1)$ecdf
    F.hat = stats::ecdf(X1)

    # Estimate the mixture model, without distinguishing each component
    mixture_hat = KDE_mixture_density(X=X2, Y=Y, null_cdf=F.hat, ker=ker)

    # Estimate the proportion of outliers in the augmented test set
    # pi.not = estimate_propOut_Storey(X=X2,Y=Y)
    # pi.not = estimate_mixing_prop(X=FX2, Y=FY, F_null=stats::punif)
    # pi.not = estimate_mixing_prop(X=X2, Y=Y, F_null=F.hat)
    pooled = c(X2,Y)
    pi.not = mixmodel::mix.model(F.hat(pooled), method = "fixed",
                        c.n = .05*log(log(length(pooled))), gridsize = 600)$alp.hat
    # Extrapolate estimate of g
    g_hat = invert_mixture(mixture_density=mixture_hat,
                           null_density=stats::dunif, prop.out=pi.not)
  } else {
    stopifnot("Error: constraint must be either increasing, decreasing"= constraint%in%c("decreasing", "increasing"))

    # pooled = c(X2,Y)
    # F_hat = stats::ecdf(X1)
    # est.fixed <- mixmodel::mix.model(F_hat(pooled), method = "fixed", c.n = .05*log(log(length(pooled))), gridsize = 600)
    # dec.dens = ifelse(constraint=="decreasing", TRUE, FALSE)
    # out = mixmodel::den.mix.model(est.fixed, dec.density = dec.dens)
    # out$x[1] <- 0
    # out$x <- out$x[-length(out$x)]
    # out$y <- out$y[-length(out$y)]
    # g_hat = function(u) sapply(u, function(i)
    #   out$y[sum( i >= out$x )] )

    pooled = c(X2,Y)
    F_hat = stats::ecdf(X1)
    est.fixed <- mixmodel::mix.model(F_hat(pooled), method = "fixed", c.n = .05*log(log(length(pooled))), gridsize = 600)
    dec.dens = ifelse(constraint=="decreasing", TRUE, FALSE)
    out = mixmodel::den.mix.model(est.fixed, dec.density = dec.dens)
    x = out$x[-length(out$x)]
    y = out$y[-length(out$y)]
    ind.min = which(x==min(x))
    ind.max = which(x==max(x))
    yleft = y[ind.min[length(ind.min)]]
    yright = y[ind.max[1]]

    g_hat = stats::approxfun(out$x[-length(out$x)],out$y[-length(out$y)], yleft=yleft, yright=yright)
  }

  return(g_hat)

}






# --------------------------------------------------------------------------------------- #
#  Implementation of the test proposed by Shiraishi (1985) with asymptotical distribution #
# --------------------------------------------------------------------------------------- #


# Stima MC delle statististiche test per singoli test points

#' stats_G_j_MC
#'
#' @param N : pooled score sample size
#' @param g : outlier distribution density
#' @param B : number of Monte Carlo repetitions used to estimate the test statistic
#'
#' @return A vector of length N corresponding to the elementary statistics of Shiraishi test (1985).
#' @export
#'
#' @examples
#' g2 = function(x, k=2) ifelse(x<1 & x>0, k*x^(k-1), 0)
#' stats_G = stats_G_j_MC(N=1000, g=g2, B=10^3)
#'
stats_G_j_MC = function(N, g, B){
  aN_j = apply(replicate(B, sapply(X=sort(stats::runif(N)), FUN=g)), 1, mean)
  return(aN_j)
}



#' stat.G
#'
#' @param Z : pooled score vector
#' @param m : calibration size
#' @param stats_G_vector : vector of Shiraishi (1985) test statistics for a fixed test sample size
#'
#' @return vector of Shiraishi test statistics given a pooled score vector
#' @export
#'
#' @examples
#' g2 = function(x, k=2) ifelse(x<1 & x>0, k*x^(k-1), 0)
#' rg2 = function(rnull, k=2) max(rnull(k))
#' stats_G = stats_G_j_MC(N=1000, g=g2, B=500)
#' stat.G = function(Z=c(runif(50),replicate(rg2(runif), 50)), m=50,stats_G_vector=stats_G)
#'
stat.G = function(Z,m,stats_G_vector){

  m = as.double(m)
  N = as.double(length(Z))

  R = rank(Z)[(m+1):N]
  T.G_i = stats_G_vector[R]

  return(T.G_i)
}


calc.stat.G = function(Z,m,stats_G_vector){

  T.G_i = stat.G(Z=Z,m=m,stats_G_vector=stats_G_vector)
  T.G_test = sum(T.G_i)

  return(T.G_test)
}


#' meanG
#'
#' @param n : test size
#' @param stats_G_vector : vector of elementary test statistics to perform the test in Shiraishi (1985). If NULL it will be computed in d_t using B_MC iterations
#'
#' @return A numeric value corresponding to the mean of the asymptotic distribution of Shiraishi test (1985)
#' @export
#'
#' @examples
#' n = 500
#' g2 = function(x, k=2) ifelse(x<1 & x>0, k*x^(k-1), 0)
#' stats_G = stats_G_j_MC(N=1000, g=g2, B=10^3)
#' mu = meanG(n=n, stats_G_vector = stats_G)
#'
meanG = function(n, stats_G_vector){

  out <- n*mean(stats_G_vector)

  return(out)
}


#' varG
#'
#' @param n : test size
#' @param m : calibration size
#' @param stats_G_vector : vector of elementary test statistics to perform the test in Shiraishi (1985). If NULL it will be computed in d_t using B_MC iterations
#'
#' @return A numeric value corresponding to the variance of the asymptotic distribution of Shiraishi test (1985)
#' @export
#'
#' @examples
#' g2 = function(x, k=2) ifelse(x<1 & x>0, k*x^(k-1), 0)
#' stats_G = stats_G_j_MC(N=1000, g=g2, B=10^3)
#' Var = varG(n=500, m=500, stats_G_vector = stats_G)
#'
varG = function(n, m, stats_G_vector){

  N = m+n
  mm <- mean(stats_G_vector)
  out = (m*n*sum((stats_G_vector - mm)^2))/(N*(N-1))
  return(out)
  #(m*h*sum((stats_G[[n-h+1]] - mm)^2))/(N*(N-1))
}




#' asymptotic.pvalue.G
#'
#' @param m : calibration size
#' @param n : test size
#' @param stats_G_vector : vector of elementary test statistics to perform the test in Shiraishi (1985). If NULL it will be computed in d_t using B_MC iterations
#' @param T.obs : observed value of the test statistic
#'
#' @return asymptotic *p*-value corresponding to the global null

asymptotic.pvalue.G <- function(m, n, stats_G_vector, T.obs) {

  mean.TG = meanG(n=n, stats_G_vector=stats_G_vector)
  variance.TG = varG(n=n, m=m, stats_G_vector=stats_G_vector)

  p.value = stats::pnorm(q=T.obs, mean=mean.TG, sd = sqrt(variance.TG), lower.tail = F)

  return(p.value)
}





#' asymptotic.critical.G
#'
#' @param m : calibration size
#' @param n : test size
#' @param stats_G_vector : vector of elementary test statistics to perform the test in Shiraishi (1985). If NULL it will be computed in d_t using B_MC iterations
#' @param alpha : significance level
#'
#' @return asymptotic critical value for the Shiraishi test statistic at level \eqn{\alpha}

asymptotic.critical.G <- function(m, n, stats_G_vector, alpha=0.1) {

  mean.TG = meanG(n=n, stats_G_vector=stats_G_vector)
  variance.TG = varG(n=n, m=m, stats_G_vector=stats_G_vector)

  critical.value = as.double(stats::qnorm(alpha, mean=mean.TG, sd = sqrt(variance.TG), lower.tail = F))

  return(critical.value)
}



#' k_mom_beta
#' @description This function computes the \eqn{k}th moment of a Beta(a,b)
#'
#' @param a : first parameter of Beta distribution
#' @param b : second parameter of Beta distribution
#' @param k : order of the moment to be computed
#'
#' @return A real number
k_mom_beta = function(a, b, k){
  den = a+b+(0:(k-1))
  num = a+(0:(k-1))
  return(prod(num/den))
}
