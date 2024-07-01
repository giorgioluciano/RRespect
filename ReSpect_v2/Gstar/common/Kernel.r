kernel <- function(H, w, s) {
  ns <- length(s)
  hs <- numeric(ns)
  
  # Uses trapezoidal rule for integration
  hs[1] <- 0.5 * log(s[2] / s[1])
  hs[ns] <- 0.5 * log(s[ns] / s[ns - 1])
  hs[2:(ns - 1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns - 2)]))
  
  S <- outer(w, s, "*")
  ws2 <- S^2
  
  K <- c((ws2 / (1 + ws2)), (S / (1 + ws2))) %*% (hs * exp(H))
  
  return(K)
}

# Example usage
H <- runif(10)  # Example substituted CRS
w <- seq(1, 10, length.out = 10)  # Example frequencies
s <- seq(1, 10, length.out = 10)  # Example relaxation modes

K <- kernel(H, w, s)
print(K)
