# Function PlotMaxwellModes
#
# Plots and compares experimental G(t) with that obtained from the
# corresponding DRS
#
# Input: g, tau = spectrum (array)
#        t  = n*1 vector contains times,
#        Gt = n*1 vector contains G(t),
#

PlotMaxwellModes <- function(g, tau, t, Gt) {
  
  N <- length(g)
  
  # Create the kernel matrix
  K <- outer(t, tau, function(T, S) exp(-T / S))
  
  # Compute the model G(t) from the Maxwell modes
  GtM <- K %*% g
  
  # Plot the experimental and model G(t)
  plot(t, Gt, log = "xy", col = "blue", pch = 16, xlab = "Time", ylab = "G(t)", main = "Comparison of Experimental and Model G(t)")
  lines(t, GtM, col = "black", lwd = 2)
}

# Example usage
# Assuming g, tau, t, and Gt are defined elsewhere in your R script
# g <- ... # Define the spectrum array
# tau <- ... # Define the relaxation times array
# t <- ... # Define the times array
# Gt <- ... # Define the experimental G(t) array
# PlotMaxwellModes(g, tau, t, Gt)
