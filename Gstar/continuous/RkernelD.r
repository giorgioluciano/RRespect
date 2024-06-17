# Function: kernelD
# Output: DK = Jacobian of H
# Input: H = substituted CRS,
#        w = n*1 vector contains frequencies,
#        s = relaxation modes,

kernelD <- function(H, w, s) {
  ns <- length(s)
  hs <- numeric(ns)
  hs[1] <- 0.5 * log(s[2] / s[1])
  hs[ns] <- 0.5 * log(s[ns] / s[ns - 1])
  hs[2:(ns - 1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns - 2)]))
  
  n <- length(w)
  
  # Create matrices S and W
  S <- matrix(rep(s, each = n), nrow = ns, byrow = TRUE)
  W <- matrix(rep(w, times = ns), nrow = ns, byrow = TRUE)
  
  ws <- S * W
  ws2 <- ws^2
  
  # Calculate Hsuper
  Hsuper <- matrix(rep(hs * exp(H), each = 2*n), nrow = 2*n, byrow = TRUE)
  
  # Calculate DK
  DK <- cbind((ws2 / (1 + ws2)), (ws / (1 + ws2))) * Hsuper
  
  return(DK)
}

# Example usage
H <- c(rep(0, length(s)))  # Example H vector (all zeros)
w <- seq(1, 10, length.out = 20)  # Example frequencies
s <- seq(0.1, 1, length.out = 10)  # Example relaxation modes

DK <- kernelD(H, w, s)
print(DK)

