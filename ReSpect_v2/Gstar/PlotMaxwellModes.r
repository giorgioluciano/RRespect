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
  X <- matrix(rep(t, each = length(w)), nrow = length(w), ncol = N)
  Y <- matrix(rep(w, N), nrow = length(w), ncol = N)
  ws <- X * Y
  ws2 <- ws^2

  # Compute GpM and GppM
  
  
  GpM_nog <- t(ws2/(1+ws2))
  GpM <- crossprod(GpM_nog,g)
  
  GppM_nog <- t(ws/(1+ws2))
  GppM <- crossprod(GppM_nog,g)
  
  # Plotting using log-log scale
  plot(w, Gp, log = "xy", col = "blue", pch = 1, xlab = "Frequency", ylab = "Modulus", type = "b", lwd = 2)
  lines(w, GpM, col = "black", lwd = 2, type = "b")
  points(w, Gpp, col = "blue", pch = 1, lwd = 2)
  lines(w, GppM, col = "black", lwd = 2, type = "b")
}

# Example usage
# Define example input values
#g <- c(1, 2, 3)
#t <- c(0.1, 0.2, 0.3)
#w <- seq(0.1, 10, length.out = 100)
#Gp <- runif(100, min = 0, max = 10)
#Gpp <- runif(100, min = 0, max = 10)

# Call the function
#PlotMaxwellModes(g, t, w, Gp, Gpp)
