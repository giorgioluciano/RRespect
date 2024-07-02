# Function: kernel
#
# Outputs the 2n*1 dimensional vector K(H)(w) which is comparable to Gexp
# Modifying kernel for unevenly spaced s_i
#
# Input: H = substituted CRS,
#        w = n*1 vector containing frequencies,
#        s = relaxation modes
#

kernel <- function(H, w, s) {
  ns <- length(s)
  hs <- numeric(ns)
  
  # Uses trapezoidal rule for integration
  hs[1] <- 0.5 * log(s[2] / s[1])
  hs[ns] <- 0.5 * log(s[ns] / s[ns - 1])
  
  hs[2:(ns - 1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns - 2)]))
  
  # Create meshgrid equivalent
  ws <- outer(w, s, "*")
  ws2 <- ws^2
  
  # Calculate K
  K <- c((ws2 / (1 + ws2)) %*% (hs * exp(H)), (ws / (1 + ws2)) %*% (hs * exp(H)))
  
  return(K)
}
