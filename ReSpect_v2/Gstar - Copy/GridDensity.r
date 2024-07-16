# Function: GridDensity(input)
#
# Takes in a PDF or density function, and spits out a bunch of points in
# accordance with the PDF
#
# Input:
#  x  = vector of points. It need *not* be equispaced,
#  px = vector of same size as x: probability distribution or
#       density function. It need not be normalized but has to be positive.
#  N  = Number of points >= 3. The end points of "x" are included
#       necessarily,
#  Pt = Optional argument. If present, then some plotting.
# 
# Output:
#  z  = Points distributed according to the density
#  hz = width of the "intervals" - useful to apportion domain to points
#       if you are doing quadrature with the results, for example.
# (c) Sachin Shanbhag, March 5, 2012
library(splines)

GridDensity <- function(x, px, N, Pt = FALSE) {
 
  #

  npts <- 100
  xi <- seq(min(x), max(x), length.out = npts)
  pint <- spline(x, px, xout = xi)$y
  # pint = interp1(x,px,xi,'spline')
  
  
  ci <- cumtrapz(xi,pint)
  pint <- pint / ci[npts]
  ci <- ci /ci[npts]

  alfa <- 1 / (N - 1)
  zij <- numeric(N)
  z <- numeric(N)

  z[1] <- min(x)
  z[N] <- max(x)

  beta <- seq(0.5, N - 1.5) * alfa
  zij <- c(z[1], approx(ci, xi, beta)$y, z[N])
  h <- diff(zij)

  beta <- seq(1, N - 2) * alfa
  z[2:(N - 1)] <- approx(ci, xi, beta)$y
  
  

  if (Pt) {
    par(mfrow = c(2, 1))

    plot(xi, ci, type = "l", col = "blue", lwd = 2, xlab = "x", ylab = "CDF(x)/PDF(x)", main = "Visualization")
    lines(xi, pint, col = "red", lwd = 2)
    legend("topright", legend = c("CDF", "PDF"), col = c("blue", "red"), lwd = 2)
	windows() # change depending on your OS
    plot(z, rep(0.5, length(z)), col = "red", pch = 19, xlab = "z", ylab = "", main = "Points Distributed According to Density")
    for (i in 2:N) {
      lines(c(zij[i], zij[i]), c(0, 1), col = "blue")
    }
  }

  return(list(z = z, hz = h))
}

# Example usage
#x <- seq(0, 10, length.out = 100)
#px <- dnorm(x, mean = 5, sd = 1)
#N <- 10
#Pt <- TRUE

#result <- GridDensity(x, px, N, Pt)
#print(result)
