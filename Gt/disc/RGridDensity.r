# Function: GridDensity
# Takes in a PDF or density function, and generates points in accordance with the PDF
# Input:
#   x  = vector of points. It need *not* be equispaced,
#   px = vector of same size as x: probability distribution or density function. It need not be normalized but has to be positive.
#   N  = Number of points >= 3. The end points of "x" are included necessarily,
#   Pt = Optional argument. If present, then some plotting.
# Output:
#   z  = Points distributed according to the density
#   h  = width of the "intervals" - useful to apportion domain to points if you are doing quadrature with the results, for example.
# (c) Sachin Shanbhag, November 11, 2015

GridDensity <- function(x, px, N, Pt = FALSE) {
  
  npts <- 100                           # can potentially change
  xi <- seq(min(x), max(x), length.out = npts)   # reinterpolate on equi-spaced axis
  pint <- spline(x, px, xi)$y          # smoothen using splines
  ci <- cumtrapz(xi, pint)             # cumulative integral
  pint <- pint / ci[length(ci)]        # normalize pint
  ci <- ci / ci[length(ci)]            # normalize ci
  
  alfa <- 1 / (N - 1)                  # alfa/2 + (N-1)*alfa + alfa/2
  zij <- rep(0, N)                     # quadrature interval end marker
  z <- rep(0, N)                       # quadrature point
  
  z[1] <- min(x)
  z[N] <- max(x)
  
  # ci(Z_j,j+1) = (j - 0.5) * alfa
  beta <- seq(0.5, N - 1.5) * alfa
  zij <- c(z[1], spline(ci, xi, beta)$y, z[N])
  h <- diff(zij)
  
  # Quadrature points are not the centroids, but rather the center of masses of the quadrature intervals
  beta <- seq(1, N - 2) * alfa
  z[2:(N-1)] <- spline(ci, xi, beta)$y
  
  # Plotting if required
  if (Pt) {
    par(mfrow = c(2, 1))
    
    plot(xi, ci, type = 'l', col = 'blue', lwd = 2, ylim = c(min(ci), max(ci)))
    lines(xi, pint, col = 'red', lwd = 2)
    legend('topright', legend = c('CDF', 'PDF'), col = c('blue', 'red'), lty = 1, lwd = 2)
    title(main = 'Visualization')
    xlabel <- 'x'
    ylabel <- 'CDF(x)/PDF(x)'
    grid()
    
    plot(z, rep(0.5, length(z)), type = 'o', col = 'red', ylim = c(0, 1), xlab = 'z')
    for (i in 2:N) {
      lines(c(zij[i], zij[i]), c(0, 1), col = 'blue')
    }
    grid()
  }
  
  return(list(z = z, h = h))
}

# Example usage:
# Replace 'x', 'px', 'N', and optionally 'Pt' with appropriate vectors/values for your data
x <- seq(0, 10, length.out = 100)  # Example vector of points
px <- dnorm(x, mean = 5, sd = 1)   # Example probability distribution function
N <- 10                            # Example number of points

# Call the GridDensity function
result <- GridDensity(x, px, N, Pt = TRUE)
print(result$z)  # Display the points distributed according to the density


  
  
