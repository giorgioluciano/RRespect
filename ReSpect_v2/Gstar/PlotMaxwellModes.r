#
# Function PlotMaxwellModes(input)
#
# Plots and compares experimental G(t) with that obtained from the
# corresponding corresponding DRS
#
# Input: g, tau = spectrum (array)
#        t  = n*1 vector contains times,
#        Gt = n*1 vector contains G(t),
#


PlotMaxwellModes <- function(g, tau, t, Gt) {
  N <- length(g)
  
  S <- outer(t, tau, FUN = function(x, y) x)
  T <- outer(t, tau, FUN = function(x, y) y)
  
  K <- exp(-T / S)
  GtM <- K %*% g
  
  plot(t, Gt, log = "xy", col = "blue", pch = 16, lwd = 2, xlab = "t", ylab = "G(t)")
  lines(t, GtM, col = "black", lwd = 2)
}
