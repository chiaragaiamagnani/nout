# Robustness analysis (RA)



RA.generate_scores = function(m, n, theta, gamma, k, nu, F.norm = T){

  # Vector that indicates if a test observation is an inlier (0) or an outlier (1)
  b_theta = rbinom(n=n, size=1, prob=theta)

  n.inlier = sum(b_theta); n.outlier = n-n.inlier
  S_X = rnorm(n=n.inlier)

  # Vector that indicates if an outlier observation in the test set
  # is drawn from F^k (1) or from N(mu.F-nu, var.F) (0), where mu.F=0 and var.F=1
  b_gamma = rbinom(n=n.outlier, size=1, prob=gamma)

  n1.Fk = sum(b_gamma); n1.Nshifted = n1-n1.Fk

  S_Y.Nshifted = rnorm(n=n1.Nshifted, mean=nu)
  S_Y.Fk = rep(max(rnorm(n=k)), times=n1.Fk)
  S_Y = c(S_Y.Fk, S_Y.Nshifted)


  # In this case F is not a standard Normal N(0,1), but a standard Uniform U(0,1)
  if(!F.norm){

    # Vector that indicates if a test observation is an inlier (0) or an outlier (1)
    b_theta = rbinom(n=n, size=1, prob=theta)

    n.inlier = sum(b_theta); n.outlier = n-n.inlier
    S_X = runif(n=n.inlier)

    # Vector that indicates if an outlier observation in the test set
    # is drawn from F^k (1) or from N(mu.F-nu, var.F) (0), where mu.F=0 and var.F=1
    b_gamma = rbinom(n=n.outlier, size=1, prob=gamma)

    n1.Fk = sum(b_gamma); n1.Nshifted = n1-n1.Fk

    S_Y.Nshifted = rnorm(n=n1.Nshifted, mean=nu)
    S_Y.Fk = rep(max(runif(n=k)), times=n1.Fk)
    S_Y = c(S_Y.Fk, S_Y.Nshifted)

  }

  return(list("S_X" = S_X,
              "S_Y" = S_Y,
              "m" = m,
              "n" = n,
              "theta" = theta,
              "gamma" = gamma,
              "k" = k,
              "nu" = nu,
              "F.norm" = F.norm))

}




RA.compare_methods = function(S_X, S_Y, alpha){

  d_BH = nout::d_benjhoch(S_Y=S_Y, S_X=S_X, alpha=alpha)
  d_T2 = nout::d_MannWhitney(S_Y=S_Y, S_X=S_X, alpha=alpha)
  d_T3 = nout::d_MannWhitneyk3(S_Y=S_Y, S_X=S_X, alpha=alpha)
  d_S = nout::d_Simes(S_Y=S_Y, S_X=S_X, alpha=alpha)
  d_F = nout::d_Fisher.corrected(S_Y=S_Y, S_X=S_X, alpha=alpha)

  return(list("d_BH"=d_BH,
              "d_T2"=d_T2,
              "d_T3"=d_T3,
              "d_S"=d_S,
              "d_F"=d_F))

}



