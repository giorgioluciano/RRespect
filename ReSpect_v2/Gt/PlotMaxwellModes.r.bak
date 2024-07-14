PlotMaxwellModes <- function(g, tau, t, Gt) {
  # Function PlotMaxwellModes(input)
  #
  # Plots experimental dynamic moduli and dynamic moduli obtained from the
  # corresponding CRS
  #
  # Input: g, t = spectrum
  #        w = n*1 vector contains frequencies,
  #        Gp, Gpp
  #

  N <- length(g)

  # Create a grid of t and w
  res <- meshgrid(tau,t)
  
  S =res$X
  T= res$Y
  
  
  K <- (exp(-T/S))
  GtM <- K %*% g
  plot(t, Gt, log = "xy", pch = 16, col = "blue", xlab = "t", ylab = "Gt", main = "Log-log plot")
  lines(t, GtM, col = "black", lwd = 2)
  
  # Aggiunta della legenda
  legend("topright", legend = c("Gt", "GtM"), col = c("blue", "black"), lty = c(NA, 1), pch = c(16, NA), lwd = c(NA, 2))
}