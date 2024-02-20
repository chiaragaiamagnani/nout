

#' d_selection_t
#'
#' @description It returns the lower bound for the number of true discoveries in an arbitrary selection set \eqn{S}, which is a
#' subset of the test set, obtained with closed testing procedure using LMPI \eqn{T_2} (Wilcoxon-Mann-Whitney)
#' or LMPI \eqn{T_3} or Fisher local test.
#'
#' @param S_Y : test score vector
#' @param S_X : calibration score vector
#' @param S : set of selected indices in the test set on which compute the lower bound for the number of outliers
#' @param statistic : parameter indicating the local test to be used in closed testing procedure.
#' It can be either \eqn{T_2, T_3} or adjusted Fisher test
#' @param alpha : significance level
#' @param n_perm : if \eqn{min(m,n)\leq n_perm} critical values will be computed via permutation. Default value is 10
#' @param B : number of permutation to compute critical values. Default value is 10^3
#' @param critical_values : if not \code{NULL}, a vector of precomputed critical values obtained using
#' the permutation distribution of the test statistic
#' @param seed : seed to ensure reproducible results
#'
#' @return A list:
#' \itemize{
#' \item \code{lower_bound}: an integer which is the \eqn{(1 − \alpha)}-confidence lower bound for
#' the number of true discoveries in the selected set \eqn{S} in closed testing procedure using the chosen local test
#' \item \code{S}: the selection set, i.e., the selected subset of the test indices
#' \item \code{global.pvalue}: the global *p*-value, i.e., the *p*-value that closed testing procedure uses to reject the global null
#' }
#'
#' @export
#'
#' @examples
#' set.seed(321)
#' Sxy = sample(x=1:1000, size=100)
#' Sx = sample(Sxy, size=70); m=length(Sx)
#' Sy = setdiff(Sxy, Sx); n=length(Sy)
#' d_t(S_Y=Sy, S_X=Sx, statistic="T2", alpha=0.1)
#' d_selection_t(S_Y=Sy, S_X=Sx, statistic="T2", alpha=0.1)
#' d_t(S_Y=Sy, S_X=Sx, statistic="T3", alpha=0.1)
#' d_selection_t(S_Y=Sy, S_X=Sx, statistic="T3", alpha=0.1)
#' d_t(S_Y=Sy, S_X=Sx, statistic="fisher", alpha=0.1)
#' d_selection_t(S_Y=Sy, S_X=Sx, statistic="fisher", alpha=0.1)
#'
#' d_selection_t(S_Y=Sy, S_X=Sx, S=c(1:3,17), statistic="T2", alpha=0.1)
#' d_selection_t(S_Y=Sy, S_X=Sx, S=c(1:3,17, 7:11, 23, 21), statistic="T2", n_perm=8, alpha=0.1)
#' d_selection_t(S_Y=Sy, S_X=Sx, S=c(1,3,7,19), statistic="T3", alpha=0.1)
#' d_selection_t(S_Y=Sy, S_X=Sx, S=c(1,3,7,19), statistic="fisher", alpha=0.1)
#' d_selection_t(S_Y=Sy, S_X=Sx, S=c(1,3,7,23,11,28,19), n_perm=6, statistic="fisher", alpha=0.1)

d_selection_t <- function(S_Y, S_X, S=NULL, statistic="T2", alpha=0.1, n_perm=10, B=10^3, critical_values=NULL, seed=123){

  stopifnot(statistic %in% c("T2", "T3", "fisher"))

  if(statistic=="T2") {
    stat.func = stat.T2
    asymptotic.critical.func = asymptotic.critical.T2
    asymptotic.pvalue.func = asymptotic.pvalue.T2
    } else if (statistic=="T3") {
    stat.func = stat.T3
    asymptotic.critical.func = asymptotic.critical.T3
    asymptotic.pvalue.func = asymptotic.pvalue.T3
  } else if (statistic=="fisher") {
    stat.func = stat.Fisher
    asymptotic.critical.func = asymptotic.critical.Fisher
    asymptotic.pvalue.func = asymptotic.pvalue.Fisher
  }

  n = length(S_Y)
  m = length(S_X)

  if(!is.null(S)){
      S = as.vector(S)
  }

  # Compute individual statistics for each test point
  S_Z = c(S_X, S_Y)
  R = stat.func(S_Z, m)

  # Compute all critical values for (m,k) from k in {1,...,n}
  crit = compute.critical.values(m, n, alpha, stat.func, asymptotic.critical.func, n_perm=n_perm, B=B, critical_values=critical_values, seed=seed)

  # Compute lower bound for S
  res = sumSome::sumStatsPar(g = R, S = S, alpha = alpha, cvs = crit)

  ## Compute p-value for the global null
  T.global = sum(R)
  pval.global = compute.global.pvalue(T.obs=T.global, m=m, n=n, stat.func=stat.func,
                                      asymptotic.pvalue.func=asymptotic.pvalue.func, n_perm=n_perm, B=B, seed=seed)

  ## Compute p-value for the selected null
  ## NOTE: this calculation is missing
  pval.selection = 1

  out = list("lower.bound" = res$TD,
             "global.pvalue" = pval.global,
             "S" = S,
             "selection.p.value" = 1)


  return(out)
}
