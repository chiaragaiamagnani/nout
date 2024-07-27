# Function to sample data from specific distributions
#
# Args:
#   n: Integer. The number of samples to generate.
#   distribution: String. The distribution to sample from.
#                 Supported values are "uniform" and "lehmann-kX" where X is a numeric value.
#
# Returns:
#   A numeric vector of sampled values.
sample_scores <- function(n, distribution) {
  if(n==0) return(NULL)
  
  # Helper function to generate uniform random numbers
  rg_null <- function(n) { runif(n) }
  
  # Initialize output vector
  Z <- numeric(n)
  
  if (distribution == "uniform") {
    # Sample from a uniform distribution
    Z <- rg_null(n)
  } else if (startsWith(distribution, "lehmann-k")) {
    # Extract the value of k from the distribution string
    k <- as.numeric(str_split(distribution, "-k")[[1]][2])
    # Sample from a Lehmann distribution with parameter k
    Z <- replicate(n, max(rg_null(k)))
  } else {
    stop("Error: unknown distribution.")
  }
  
  return(Z)
}


# Function to calculate the probability density function for specific distributions
#
# Args:
#   x: Numeric vector. The values at which to calculate the density.
#   distribution: String. The distribution for the density calculation.
#                 Supported values are "uniform" and "lehmann-kX" where X is a numeric value.
#
# Returns:
#   A numeric vector of density values corresponding to x.
density_scores <- function(x, distribution) {
  # Helper function to generate uniform random numbers
  rg_null <- function(n) { runif(n) }
  
  # Initialize output vector
  out <- numeric(length(x))
  
  if (distribution == "uniform") {
    # Calculate density for a uniform distribution
    out <- ifelse(x >= 0 & x <= 1, 1, 0)
  } else if (startsWith(distribution, "lehmann-k")) {
    # Extract the value of k from the distribution string
    k <- as.numeric(str_split(distribution, "-k")[[1]][2])
    # Calculate density for a Lehmann distribution with parameter k
    out <- ifelse(x >= 0 & x <= 1, k * x^(k - 1), 0)
  } else {
    stop("Error: unknown distribution.")
  }
  
  return(out)
}


# Function to generate a calibration set of inliers from a uniform distribution
# and a test set containing a mixture of inliers and outliers from a different distribution.
#
# Args:
#   n_cal: Integer. The number of samples in the calibration set.
#   n_test: Integer. The number of samples in the test set.
#   prop_out: Numeric. The proportion of outliers in the test set.
#   alternative: String. The type of distribution for outliers in the test set.
#
# Returns:
#   A list containing:
#     - scores.cal: Numeric vector. Calibration set scores from a uniform distribution.
#     - scores.test: Numeric vector. Test set scores with a mixture of inliers and outliers.
#     - outlier.test: Integer vector. Indicator of outliers in the test set (0 = inlier, 1 = outlier).
generate_cal_test_scores <- function(n_cal, n_test, prop_out, alternative) {
  # Generate calibration set scores from a uniform distribution
  Z_cal <- sample_scores(n_cal, "uniform")
  
  # Calculate the number of outliers and inliers in the test set
  n_test_1 <- as.integer(round(prop_out * n_test))
  n_test_0 <- n_test - n_test_1
  
  # Generate inlier test set scores from a uniform distribution
  if(n_test_0 > 0) {
    Z_test_0 <- sample_scores(n_test_0, "uniform")
  } else {
    
  }
  
  # Generate outlier test set scores from the specified alternative distribution
  Z_test_1 <- sample_scores(n_test_1, alternative)
  
  # Combine inlier and outlier scores to create the test set
  Z_test <- c(Z_test_0, Z_test_1)
  
  # Create an indicator vector for outliers in the test set
  outlier.test <- c(rep(0, n_test_0), rep(1, n_test_1))
  
  # Randomly shuffle the order of the test set scores and outlier indicators
  order <- sample(n_test)
  Z_test <- Z_test[order]
  outlier.test <- outlier.test[order]
  
  # Return the calibration set scores, test set scores, and outlier indicators as a list
  out <- list(scores.cal = Z_cal, scores.test = Z_test, outlier.test = outlier.test)
  return(out)
}