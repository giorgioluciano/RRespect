library(minpack.lm)

# Function to find parameters with nlsLM and optional exclusion of outliers or points from head/tail
find_parameters_lm <- function(file_path, exclude_head = 0, exclude_tail = 0) {
  # Step 1: Read the data from the file
  data <- read.table(file_path, header = FALSE)
  gamma <- data$V1  # Extract the gamma (shear stress) values
  h <- exp(data$V2)      # Extract the h (damping function) values
  
  # Step 2: Manually exclude points from the head (start) and tail (end)
  if (exclude_head > 0) {
    gamma <- gamma[-(1:exclude_head)]  # Remove the first 'exclude_head' points
    h <- h[-(1:exclude_head)]
  }
  if (exclude_tail > 0) {
    gamma <- gamma[-((length(gamma) - exclude_tail + 1):length(gamma))]  # Remove the last 'exclude_tail' points
    h <- h[-((length(h) - exclude_tail + 1):length(h))]
  }
  
  # Step 3: Fit the model using nlsLM with bounds to ensure positive parameters
  fit <- nlsLM(h ~ 1 / (1 + a * gamma^b),
               start = list(a = 1, b = 0.5),  # Initial guesses for a and b
               #lower = c(a = 0, b = 0),     # Lower bounds: a and b must be positive
               data = data.frame(gamma, h)) # Provide the data for the model
  
  # Return the fitted model and data for further diagnostics
  return(list(fit = fit, gamma = gamma, h = h))
}

# Diagnostic function to plot log10-log10 of gamma vs h, with observed and predicted values
diagnostics_plot <- function(fit, gamma, h) {
  # Extract fitted parameters from the model
  coef_fit <- coef(fit)
  a <- coef_fit["a"]
  b <- coef_fit["b"]
  
  # Compute predicted h values based on the fitted model
  h_pred <- 1 / (1 + a * gamma^b)
  
  # Create the log-log plot of observed gamma vs h
  plot(log10(gamma), log10(h), col = "blue", pch = 19,
       xlab = expression(log[10](gamma)),  # Label for the x-axis
       ylab = expression(log[10](h)),      # Label for the y-axis
       main = "Diagnostic Plot: log10(gamma) vs log10(h)")
  
  # Add predicted values to the plot in red
  points(log10(gamma), log10(h_pred), col = "red", pch = 4)
  
  # Add a legend to distinguish observed vs predicted points
  legend("topright", legend = c("Observed", "Predicted"), col = c("blue", "red"), pch = c(19, 4))
}

# Example usage:
# - Exclude the first 2 points (head) and the last 3 points (tail)
result <- find_parameters_lm("./output/H.dat", exclude_head = 12, exclude_tail = 13)

# Generate the diagnostic plot to compare observed vs predicted values
diagnostics_plot(result$fit, result$gamma, result$h)
