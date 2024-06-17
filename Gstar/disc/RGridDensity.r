# Load necessary libraries
library(pracma)  # For the 'cumtrapz' function
library(splines) # For spline interpolation

# Function to generate points distributed according to the given density function
# Input: x  = vector of points (not necessarily equispaced)
#        px = probability distribution or density function of the same size as x
#        N  = number of points >= 3 (including end points of x)
#        plot = logical, if TRUE, will plot the results
# Output: z  = points distributed according to the density
#         hz = width of the intervals

GridDensity <- function(x, px, N, plot = FALSE) {
  npts <- 100  # Can be changed
  xi <- seq(min(x), max(x), length.out = npts)  # Reinterpolate on equi-spaced axis
  pint <- splinefun(x, px, method = "fmm")(xi)  # Smoothen using splines
  ci <- cumtrapz(xi, pint)  # Cumulative integral
  pint <- pint / ci[npts]
  ci <- ci / ci[npts]  # Normalize ci
  
  alfa <- 1 / (N - 1)  # alfa/2 + (N-1)*alfa + alfa/2
  zij <- numeric(N)  # Quadrature interval end marker
  z <- numeric(N)  # Quadrature point
  
  z[1] <- min(x)
  z[N] <- max(x)
  
  # ci(Z_j, j+1) = (j - 0.5) * alfa
  beta <- seq(0.5, N - 1.5) * alfa
  zij <- c(z[1], splinefun(ci, xi, method = "fmm")(beta), z[N])
  hz <- diff(zij)
  
  # Quadrature points are the center of masses of the quadrature intervals
  beta <- seq(1, N - 2) * alfa
  z[2:(N - 1)] <- splinefun(ci, xi, method = "fmm")(beta)
  
  # Optional plotting
  if (plot) {
    par(mfrow = c(2, 1))
    
    plot(xi, ci, type = "l", col = "blue", lwd = 2, main = "Visualization", ylab = "CDF(x)/PDF(x)", xlab = "x")
    lines(xi, pint, col = "red", lwd = 2)
    legend("topright", legend = c("CDF", "PDF"), col = c("blue", "red"), lty = 1, lwd = 2)
    
    plot(z, rep(0.5, length(z)), type = "p", col = "red", main = "Generated Points", ylab = "", xlab = "z", ylim = c(0, 1))
    abline(h = 0.5, col = "gray", lty = 2)
    for (i in 2:N) {
      lines(c(zij[i], zij[i]), c(0, 1), col = "blue")
    }
  }
  
  return(list(z = z, hz = hz))
}

# Example usage
x <- seq(0, 10, length.out = 50)
px <- dnorm(x, mean = 5, sd = 2)  # Example density function
N <- 10

result <- GridDensity(x, px, N, plot = TRUE)
print(result$z)
print(result$hz)

  
