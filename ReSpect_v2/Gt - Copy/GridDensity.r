# Function: GridDensity
#
# Takes in a PDF or density function, and spits out a bunch of points in
# accordance with the PDF
#
# Input:
#   x  = vector of points. It need *not* be equispaced,
#   px = vector of same size as x: probability distribution or
#        density function. It need not be normalized but has to be positive.
#   N  = Number of points >= 3. The end points of "x" are included necessarily,
#   Pt = Optional argument. If present, then some plotting.
# 
# Output:
#   z  = Points distributed according to the density
#   hz = width of the "intervals" - useful to apportion domain to points
#        if you are doing quadrature with the results, for example.
#
# (c) Sachin Shanbhag, November 11, 2015
library(pracma)


GridDensity <- function(x, px, N, Pt = FALSE) {

  npts <- 100                               # can potentially change
  xi <- seq(min(x), max(x), length.out = npts)  # reinterpolate on equi-spaced axis
  pint <- splinefun(x, px)(xi)              # smoothen using splines
  ci <- cumtrapz(xi, pint)
  pint <- pint / ci[length(ci)]
  ci <- ci / ci[length(ci)]                 # normalize ci
  
  alfa <- 1 / (N - 1)                       # alfa/2 + (N-1)*alfa + alfa/2
  zij <- numeric(N)                         # quadrature interval end marker
  z <- numeric(N)                           # quadrature point
  
  z[1] <- min(x)
  z[N] <- max(x)
  
  # ci(Z_j,j+1) = (j - 0.5) * alfa
  beta <- seq(0.5, N - 1.5) * alfa
  zij <- c(z[1], approx(ci, xi, xout = beta, method = "linear")$y, z[N])
  hz <- diff(zij)
  
  # Quadrature points are not the centroids, but rather the center of masses
  # of the quadrature intervals
  beta <- seq(1, N - 2) * alfa
  z[2:(N - 1)] <- approx(ci, xi, xout = beta, method = "linear")$y
  
  # Some plotting if required
  if (Pt) {
    par(mfrow = c(2, 1))
    
    plot(xi, ci, type = "l", col = "blue", lwd = 2, ylim = c(0, 1),
         xlab = "x", ylab = "CDF(x)/PDF(x)", main = "Visualization")
    lines(xi, pint, col = "red", lwd = 2)
    legend("topright", legend = c("CDF", "PDF"), col = c("blue", "red"), lwd = 2)
    
    plot(z, rep(0.5, length(z)), type = "p", col = "red", pch = 16,
         xlab = "z", ylab = "", main = "Quadrature Points")
    for (i in 2:N) {
      lines(c(zij[i], zij[i]), c(0, 1), col = "blue")
    }
  }
  
  return(list(z = z, hz = hz))
}

  
  
