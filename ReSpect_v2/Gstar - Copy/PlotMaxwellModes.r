PlotMaxwellModes <- function(g, t, w, Gp, Gpp) {
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
  res <- pracma::meshgrid(t,w)
  X <- res$X
  Y <- res$Y
  ws <- X * Y
  ws2 <- ws^2
  
  # Compute GpM and GppM
  GpM <- rowSums((ws2 / (1 + ws2)) * rep(g, each = 100))
  GppM <- rowSums((ws / (1 + ws2)) * rep(g, each = 100))
  
  # Plotting using log-log scale
  plot(w, Gp, log = "xy", col = "blue", pch = 1, xlab = "Frequency", ylab = "Modulus", type = "b", lwd = 2)
  lines(w, GpM, col = "black", lwd = 2, type = "b")
  points(w, Gpp, col = "blue", pch = 1, lwd = 2)
  lines(w, GppM, col = "black", lwd = 2, type = "b")
}