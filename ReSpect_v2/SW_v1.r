library(minpack.lm)

# Function for fitting the Soskey model
soskey_fit <- function(w, h, plot_results = TRUE) {
  
  # Define the Soskey model
  soskey_model <- function(w, a, b) {
    1 / (1 + (a * abs(w)^b))
  }
  
  # Create a data frame
  data <- data.frame(w = w, h = h)
  
  # Fit the model using non-linear least squares
  fit <- nlsLM(h ~ soskey_model(w, a, b),
               data = data,
               start = list(a = 0.9, b = 0.5))
  
  # Extract the fitted parameters
  a_fitted <- coef(fit)["a"]
  b_fitted <- coef(fit)["b"]
  
  # Print the fitted parameters
  cat("Fitted parameters:\n")
  cat("a =", a_fitted, "\n")
  cat("b =", b_fitted, "\n")
  
  # Plot the results if requested
  if (plot_results) {
    w_range <- seq(min(w), max(w), length.out = 100)
    h_fitted <- soskey_model(w_range, a_fitted, b_fitted)
    
    plot(w, h, log = "x", pch = 16, 
         xlab = "w", ylab = "h(w)", 
         main = "Soskey Model Fit")
    lines(w_range, h_fitted, col = "red")
    legend("topright", legend = c("Experimental data", "Fitted model"), 
           pch = c(16, NA), lty = c(NA, 1), col = c("black", "red"))
  }
  
  # Return the fit object for further analysis
  return(fit)
}

# # Example usage
# w <- c(0.0002710786, 0.0003097842, 0.0003525654, 0.0005018405, 0.0005584607, 
#        0.0010067603, 0.0100390038, 0.0100390038, 0.0492845550, 0.0494975295, 
#        0.1018789869, 0.2548022532, 0.2777518337, 0.5043261525, 0.7504720063, 
#        0.9769139192, 1.3481557936, 1.5144180935)
# 
# h <- c(0.95668943, 0.98463218, 0.92097764, 0.93052579, 0.88205264, 
#        0.98418234, 0.41927501, 0.36132655, 0.20219128, 0.18400445, 
#        0.14789440, 0.07450927, 0.06852493, 0.05033053, 0.03520427, 
#        0.02558643, 0.02116185, 0.0184909)
# 
# # Call the function with the plot
# fit <- soskey_fit(w, h, plot_results = TRUE)
