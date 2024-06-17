# Function: PlotMaxwellModes
# Plots and compares experimental G(t) with that obtained from the corresponding DRS
# Input:
#   g     = spectrum (array)
#   tau   = array of relaxation times
#   t     = n*1 vector contains times
#   Gt    = n*1 vector contains G(t)
PlotMaxwellModes <- function(g, tau, t, Gt) {
  
  N <- length(g)
  
  # Construct matrix K
  K <- outer(tau, t, function(tau, t) exp(-t / tau))
  
  # Compute G(t) from inferred spectrum g
  GtM <- K %*% g
  
  # Plot G(t) and G(t) from DRS
  plot(t, Gt, type = "o", col = "blue", log = "xy", xlab = "t", ylab = "G(t)", main = "Experimental vs. DRS")
  lines(t, GtM, col = "black", lwd = 2)
  
  # Add legend
  legend("topright", legend = c("Experimental G(t)", "DRS G(t)"), col = c("blue", "black"), lty = c(1, 1), lwd = c(1, 2))
}
