library(minpack.lm)

# Soskey model function
soskey_model <- function(w, a, b) {
  1 / (1 + (a * abs(w)^b))
}

# Function to fit the Soskey model with an option to plot
fit_soskey_model <- function(w, h, plot_results = TRUE) {
  
  # Prepare data for fitting
  data <- data.frame(w = w, h = h)
  
  # Fit the model using non-linear least squares
  fit <- nlsLM(h ~ soskey_model(w, a, b),
               data = data,
               start = list(a = 0.9, b = 0.5))
  
  # Extract fitted parameters
  a_fitted <- coef(fit)["a"]
  b_fitted <- coef(fit)["b"]
  
  # Generate fitted curve
  w_range <- seq(min(w), max(w), length.out = 100)
  h_fitted <- soskey_model(w_range, a_fitted, b_fitted)
  
  # Plot results if requested
  if (plot_results) {
    plot(w, h, log = "x", pch = 16,
         xlab = "w", ylab = "h(w)",
         main = "Fit of the Soskey model")
    lines(w_range, h_fitted, col = "red")
    legend("topright", legend = c("Experimental data", "Fitted model"), 
           pch = c(16, NA), lty = c(NA, 1), col = c("black", "red"))
  }
  
  # Print fitted parameters
  cat("Fitted parameters:\n")
  cat("a =", a_fitted, "\n")
  cat("b =", b_fitted, "\n")
  
  # Return the fit object for further analysis if needed
  return(list(a = a_fitted, b = b_fitted, fit = fit))
}

# Example of usage
# w <- your_w_values
# h <- your_h_values
# fit_soskey_model(w, h, plot_results = TRUE)  # Set plot_results to FALSE to exclude the plot
