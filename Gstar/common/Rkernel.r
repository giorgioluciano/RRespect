# Function: kernel
# Purpose: Computes the 2n*1 dimensional vector K(H)(w) comparable to Gexp
#          Modifies kernel for unevenly spaced s_i using trapezoidal rule for integration
# Input: H = substituted CRS,
#        w = n*1 vector containing frequencies,
#        s = relaxation modes
# Output: K = 2n*1 dimensional vector K(H)(w)
kernel <- function(H, w, s) {
  ns <- length(s)
  hs <- numeric(ns)
  
  # Calculate hs using trapezoidal rule for integration
  hs[1] <- 0.5 * log(s[2] / s[1])
  hs[ns] <- 0.5 * log(s[ns] / s[ns - 1])
  hs[2:(ns - 1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns - 2)]))
  
  # Create meshgrid for S and W
  S <- outer(s, w, "*")
  W <- outer(w, s, "*")
  
  # Compute ws and ws2
  ws <- S * W
  ws2 <- ws^2
  
  # Calculate K(H)(w)
  K <- cbind(ws2 / (1 + ws2), ws / (1 + ws2)) %*% (hs * exp(H))
  
  return(K)
}

# Example usage:
# Replace `H`, `w`, and `s` with appropriate values.
# K <- kernel(H, w, s)
# This will compute the K(H)(w) vector based on the provided H, w, and s.

