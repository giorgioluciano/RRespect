# Function: kernelD
# Output: DK = Jacobian of H
# Input: H = substituted CRS,
#        t = n*1 vector of times,
#        s = relaxation modes,
kernelD <- function(H, t, s) {
  
  ns <- length(s)
  hs <- numeric(ns)
  hs[1] <- 0.5*log(s[2]/s[1])
  hs[ns] <- 0.5*log(s[ns]/s[ns-1])
  hs[2:(ns-1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns-2)]))
  
  n <- length(t)
  
  S <- matrix(rep(s, n), nrow = n, byrow = TRUE)
  T <- matrix(rep(t, ns), nrow = n, byrow = FALSE)
  kern <- exp(-T / S)
  
  Hsuper <- t(apply(exp(H), 1, function(x) hs * x))
  DK <- kern * Hsuper
  
  return(DK)
}

# Example usage:
# Replace 't' and 's' with appropriate vectors for your data
# Replace 'H' with the actual substituted CRS vector
t <- seq(0, 10, by = 1)  # Example time vector
s <- seq(1, 5, by = 1)   # Example relaxation modes
H <- rep(0.5, length(s))  # Example substituted CRS vector

# Call the kernelD function
DK <- kernelD(H, t, s)

