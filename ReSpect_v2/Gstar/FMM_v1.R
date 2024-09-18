# Function for the generalized fractional Maxwell model with alpha and beta
fractional_maxwell <- function(omega, G0, alpha, beta) {
  complex_modulus <- G0 * (1i * omega)^alpha / (1 + (1i * omega)^beta)
  Gp <- Re(complex_modulus)  # G' (storage modulus)
  Gpp <- Im(complex_modulus) # G'' (loss modulus)
  return(list(Gp = Gp, Gpp = Gpp))
}

# Function to compute residuals between experimental and model values
residuals_fractional_maxwell <- function(params, omega, Gp_exp, Gpp_exp) {
  G0 <- params[1]
  alpha <- params[2]
  beta <- params[3]
  
  
  # Compute model G' and G''
  model <- fractional_maxwell(omega, G0, alpha, beta)
  
  # Residuals for both G' and G''
  res_Gp <- Gp_exp - model$Gp
  res_Gpp <- Gpp_exp - model$Gpp
  
  # Combine residuals into a single vector
  return(c(res_Gp, res_Gpp))
}

# Fit fractional Maxwell model to experimental data
fit_fractional_maxwell <- function(omega, Gp_exp, Gpp_exp, start_params = c(G0 = 1, alpha = 0.5, beta = 0.5)) {
  # Perform the non-linear least squares fit
  fit <- nls.lm(par = start_params, fn = residuals_fractional_maxwell, 
                omega = omega, Gp_exp = Gp_exp, Gpp_exp = Gpp_exp)
  
  # Extract fitted parameters
  fitted_params <- fit$par
  return(fitted_params)
}

# Example usage:
# Assuming `w` is the frequency, `Gp` is the storage modulus, and `Gpp` is the loss modulus
omega <- w    # Frequenza
Gp_exp <- Gp  # Modulo di conservazione sperimentale
Gpp_exp <- Gpp # Modulo di perdita sperimentale

# Fit the model to the data with initial guesses for G0, alpha, and beta
fitted_params <- fit_fractional_maxwell(omega, Gp_exp, Gpp_exp)

# Show the fitted parameters
print(fitted_params)

plot_fitting_results <- function(omega, Gp_exp, Gpp_exp, Gp_fit, Gpp_fit) {
  # Plot G' (storage modulus)
  plot(omega, Gp_exp, log = "xy", col = "blue", pch = 19, cex = 0.8,
       xlab = "Frequency (ω)", ylab = "Modulus", main = "Fractional Maxwell Model Fit")
  points(omega, Gp_fit, col = "blue", type = "l", lwd = 2)
  
  # Add G'' (loss modulus) on the same plot
  points(omega, Gpp_exp, col = "red", pch = 19, cex = 0.8)
  points(omega, Gpp_fit, col = "red", type = "l", lwd = 2)
  
  # Add legend
  legend("topright", legend = c("G' (exp)", "G' (fit)", "G'' (exp)", "G'' (fit)"),
         col = c("blue", "blue", "red", "red"), pch = c(19, NA, 19, NA), lty = c(NA, 1, NA, 1))
}

# Example usage:
# Assuming `w` is the frequency, `Gp` is the storage modulus, and `Gpp` is the loss modulus
omega <- w    # Frequenza
Gp_exp <- Gp  # Modulo di conservazione sperimentale
Gpp_exp <- Gpp # Modulo di perdita sperimentale

# Fit the model to the data with initial guesses for G0, alpha, and beta
fitted_params <- fit_fractional_maxwell(omega, Gp_exp, Gpp_exp)

# Get the fitted G' and G'' values using the optimized parameters
fitted_model <- fractional_maxwell(omega, fitted_params[1], fitted_params[2], fitted_params[3])
Gp_fit <- fitted_model$Gp
Gpp_fit <- fitted_model$Gpp

# Plot the results to see the fit
plot_fitting_results(omega, Gp_exp, Gpp_exp, Gp_fit, Gpp_fit)

# Show the fitted parameters
print(fitted_params)



