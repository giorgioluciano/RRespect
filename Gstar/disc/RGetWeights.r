# Function to find the weight of each mode by taking a weighted average of its contribution to Gp and Gpp
# Input: H = substituted CRS
#        w = n*1 vector containing frequencies
#        s = relaxation modes
# Output: wt = weight of each mode

GetWeights <- function(H, w, s) {
  ns <- length(s)
  n <- length(w)
  
  hs <- numeric(ns)
  wt <- numeric(ns)
  hs[1] <- 0.5 * log(s[2] / s[1])
  hs[ns] <- 0.5 * log(s[ns] / s[ns - 1])
  hs[2:(ns - 1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns - 2)]))
  
  S <- outer(w, s, "*")
  ws <- S
  ws2 <- ws^2
  
  wij <- rbind((ws2 / (1 + ws2)), (ws / (1 + ws2))) %*% diag(hs * exp(H))
  K <- rbind((ws2 / (1 + ws2)), (ws / (1 + ws2))) %*% (hs * exp(H))
  
  for (i in 1:(2 * n)) {
    wij[i, ] <- wij[i, ] / K[i]
  }
  
  for (j in 1:ns) {
    wt[j] <- sum(wij[, j])
  }
  
  return(wt)
}

# Example usage:
H <- rnorm(100) # Example H values
w <- seq(1, 10, length.out = 50) # Example frequency vector
s <- seq(1, 10, length.out = 100) # Example relaxation modes

weights <- GetWeights(H, w, s)
print(weights)

