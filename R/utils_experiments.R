# Function to compute global p-values using various statistical tests
#
# Args:
#   data: List. A list containing:
#     - scores.cal: Numeric vector. Calibration set scores.
#     - scores.test: Numeric vector. Test set scores.
#     - outlier.test: Numeric vector. Indicator of outliers in the test set (0 = inlier, 1 = outlier).
#   alternative: String. The type of distribution for outliers in the test set.
#
# Returns:
#   A data frame with the computed global p-values for each method.
run_global_testing <- function(data, alternative=NULL) {
  
  # Compute global p-value with Fisher's test
  pval.fisher <- d_selection_fisher(S_X = data$scores.cal, S_Y = data$scores.test, n_perm = 0, pvalue_only = TRUE)$global.pvalue
  
  # Compute global p-value with WMW test
  pval.wmw <- d_selection_higher(S_X = data$scores.cal, S_Y = data$scores.test, local.test = "WMW", n_perm = 0, pvalue_only = TRUE)$global.pvalue
  
  # Compute global p-value with WMW test (k=2)
  pval.wmw.2 <- d_selection_higher(S_X = data$scores.cal, S_Y = data$scores.test, local.test = "higher", k = 2, n_perm = 0, pvalue_only = TRUE)$global.pvalue
  
  # Compute global p-value with WMW test (k=3)
  pval.wmw.3 <- d_selection_higher(S_X = data$scores.cal, S_Y = data$scores.test, local.test = "higher", k = 3, n_perm = 0, pvalue_only = TRUE)$global.pvalue
  
  # Compute global p-value with Shirashi's approach (using oracle density)
  if(!is.null(alternative)) {
    density_oracle <- function(x) density_scores(x, alternative)
    pval.g.oracle <- d_selection_G2(S_X = data$scores.cal, S_Y = data$scores.test, g.hat = density_oracle, monotonicity = "increasing", B = 100, pvalue_only = TRUE)$global.pvalue
  } else {
    pval.g.oracle <- NA
  }
  
  # Create a data frame with the p-values
  pval_df <- tibble::tibble(
    Method = c("Fisher", "WMW", "WMW_k2", "WMW_k3", "Shirashi_oracle"),
    p.value = c(pval.fisher, pval.wmw, pval.wmw.2, pval.wmw.3, pval.g.oracle)
  )
  
  return(pval_df)
}