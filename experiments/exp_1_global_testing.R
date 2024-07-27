

# Load necessary libraries
library(tidyverse)    # For data manipulation and visualization
library(progress)     # For displaying progress bars
library(nout)         # Custom library (assuming this contains necessary functions)

# Source utility functions for data generation and experiments
source("../R/utils_data.R")
source("../R/utils_experiments.R")

###########################
## Experiment parameters ##
###########################

# Define the number of calibration data points
n_cal <- 200

# Define the number of test points
n_test <- 100

# Choose the alternative distribution for testing
#alternative <- "uniform"
alternative <- "lehmann-k2"

# List of values for the proportion of outliers in the test set
prop_out_values <- c(0, 0.25, 0.5)

# Number of repetitions for each experimental setting
n_exp <- 100

##########################
## Experiment functions ##
##########################

# Function to run a single experiment
# Args:
#   seed: Random seed for reproducibility
#   prop_out: Proportion of outliers
# Returns:
#   A tibble containing the results of the experiment
run_experiment <- function(seed, prop_out) {
  set.seed(seed)  # Set seed for reproducibility

  # Generate calibration and test data with specified parameters
  data <- generate_cal_test_scores(n_cal = n_cal, n_test = n_test, prop_out = prop_out, alternative = alternative)

  # Apply global testing methods to the generated data
  res <- run_global_testing(data, alternative = alternative)

  # Combine the results with experiment metadata (seed and proportion of outliers)
  results <- tibble(Seed = seed, Prop_Out = prop_out) %>% cbind(res)

  return(results)
}

# Function to run multiple experiments and gather results
# Args:
#   n_exp: Number of repetitions for each experimental setting
#   prop_out_values: A vector of different proportions of outliers to test
# Returns:
#   A tibble containing the combined results of all experiments
run_multiple_experiments <- function(n_exp, prop_out_values) {
  results_list <- list()  # Initialize an empty list to store results

  # Print a progress bar header
  cat("Running experiments\n")
  pb <- txtProgressBar(min = 0, max = length(prop_out_values) * n_exp, style = 3)  # Initialize progress bar

  counter <- 1  # Initialize a counter for progress tracking
  for (prop_out in prop_out_values) {  # Loop over each proportion of outliers
    for (i in 1:n_exp) {  # Loop over each repetition
      results_list[[counter]] <- run_experiment(i, prop_out)  # Run experiment and store result
      setTxtProgressBar(pb, counter)  # Update progress bar
      counter <- counter + 1  # Increment counter
    }
  }

  close(pb)  # Close the progress bar

  # Combine all results into a single tibble
  results_df <- bind_rows(results_list)

  return(results_df)
}

#####################
## Run experiments ##
#####################

# Run the experiments with specified parameters
results <- run_multiple_experiments(n_exp, prop_out_values)

##################
## Plot results ##
##################

# Define the significance level
alpha <- 0.1

# Calculate power for different methods and proportions of outliers
power_results <- results %>%
  group_by(Method, Prop_Out) %>%
  summarize(
    Power = mean(p.value < alpha),  # Calculate power as the proportion of p-values below 0.05
    SE = sqrt((Power * (1 - Power)) / n())  # Calculate standard error of the power estimate
  )

# Plot the power for different methods and proportions of outliers
power_results %>%
  ggplot(aes(x = Prop_Out, y = Power, color = Method, shape = Method)) +  # Set aesthetics
  geom_line() +  # Add lines
  geom_point() +  # Add points
  geom_errorbar(aes(ymin = Power - 2 * SE, ymax = Power + 2 * SE), width = 0.02) +  # Add error bars
  geom_hline(yintercept = alpha, linetype = 2) +  # Add horizontal line at significance level
  theme_minimal(base_size = 15) +  # Set minimal theme for the plot
  labs(
    title = "Power of Different Methods for Various Proportions of Outliers",
    subtitle = sprintf("N-cal: %d, N-test: %d, Alternative distribution: %s", n_cal, n_test, alternative),
    x = "Proportion of Outliers",
    y = "Power",
    color = "Method"
  ) +
  ylim(0, 1)  # Set y-axis limits
