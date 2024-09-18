library(minpack.lm)
library(ggplot2)

# Function for the generalized fractional Maxwell model
fractional_maxwell <- function(omega, G0, alpha, beta) {
  complex_modulus <- G0 * (1i * omega)^alpha / (1 + (1i * omega)^beta)
  Gp <- Re(complex_modulus)
  Gpp <- Im(complex_modulus)
  return(list(Gp = Gp, Gpp = Gpp))
}

# Function to compute residuals
residuals_fractional_maxwell <- function(params, omega, Gp_exp, Gpp_exp) {
  G0 <- params[1]
  alpha <- params[2]
  beta <- params[3]
  
  model <- fractional_maxwell(omega, G0, alpha, beta)
  
  res_Gp <- Gp_exp - model$Gp
  res_Gpp <- Gpp_exp - model$Gpp
  
  return(c(res_Gp, res_Gpp))
}

# Function to fit the fractional Maxwell model
fit_fractional_maxwell <- function(omega, Gp_exp, Gpp_exp, start_params = c(G0 = 1e5, alpha = 0.5, beta = 0.5)) {
  fit <- nls.lm(par = start_params, 
                fn = residuals_fractional_maxwell, 
                omega = omega, Gp_exp = Gp_exp, Gpp_exp = Gpp_exp,
                lower = c(0, 0, 0),
                upper = c(Inf, 1, 1))
  
  return(fit$par)
}

# Function to plot results using ggplot2
plot_fitting_results <- function(omega, Gp_exp, Gpp_exp, Gp_fit, Gpp_fit) {
  data <- data.frame(
    omega = rep(omega, 4),
    Modulus = c(Gp_exp, Gpp_exp, Gp_fit, Gpp_fit),
    Type = rep(c("G' (exp)", "G'' (exp)", "G' (fit)", "G'' (fit)"), each = length(omega))
  )
  
  ggplot(data, aes(x = omega, y = Modulus, color = Type)) +
    geom_point(data = subset(data, Type %in% c("G' (exp)", "G'' (exp)"))) +
    geom_line(data = subset(data, Type %in% c("G' (fit)", "G'' (fit)"))) +
    scale_x_log10() + scale_y_log10() +
    labs(x = "Frequency (ω)", y = "Modulus", title = "Fractional Maxwell Model Fit") +
    theme_minimal()
}

# Main execution
# Assuming 'dati' dataframe exists with columns 'w', 'Gp', and 'Gpp'
omega <- dati$w
Gp_exp <- dati$Gp
Gpp_exp <- dati$Gpp

# Fit the model
fitted_params <- fit_fractional_maxwell(omega, Gp_exp, Gpp_exp)

# Calculate fitted values
fitted_model <- fractional_maxwell(omega, fitted_params[1], fitted_params[2], fitted_params[3])
Gp_fit <- fitted_model$Gp
Gpp_fit <- fitted_model$Gpp

# Plot results
p <- plot_fitting_results(omega, Gp_exp, Gpp_exp, Gp_fit, Gpp_fit)
print(p)

# Print fitted parameters
cat("Fitted parameters:\n")
cat("G0:", fitted_params[1], "\n")
cat("alpha:", fitted_params[2], "\n")
cat("beta:", fitted_params[3], "\n")
