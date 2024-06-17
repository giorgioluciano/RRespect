# Function: GetWeights
# Finds the weight of "each" mode by taking a weighted average of its contribution
# to G(t)
# Input: H = CRS,
#        t = n*1 vector contains times,
#        s = relaxation modes,
# Output: wt = weight of each mode
GetWeights <- function(H, t, s) {
  
  ns <- length(s)
  n <- length(t)
  
  hs <- rep(0, ns)
  wt <- rep(0, ns)
  
  hs[1] <- 0.5 * log(s[2] / s[1])
  hs[ns] <- 0.5 * log(s[ns] / s[ns - 1])
  hs[2:(ns - 1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns - 2)]))
  
  # Compute kernel matrix
  S <- matrix(rep(s, n), n, ns, byrow = TRUE)
  T <- matrix(rep(t, each = ns), n, ns, byrow = TRUE)
  kern <- exp(-T / S)
  
  # Compute wij
  wij <- kern * diag(hs * exp(H))
  
  # Compute K
  K <- kern %*% (hs * exp(H))
  
  # Normalize wij
  for (i in 1:n) {
    wij[i, ] <- wij[i, ] / K[i]
  }
  
  # Compute weights wt
  for (j in 1:ns) {
    wt[j] <- sum(wij[, j])
  }
  
  return(wt)
}

# Example usage:
# Replace 'H', 't', and 's' with appropriate vectors/values for your data
H <- runif(5)     # Example CRS vector
t <- seq(0, 10, length.out = 100)  # Example time vector
s <- seq(1, 5, by = 1)             # Example relaxation modes

# Call the GetWeights function
weights <- GetWeights(H, t, s)
print(weights)

