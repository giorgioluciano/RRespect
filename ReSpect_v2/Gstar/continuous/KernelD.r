kernelD <- function(H, w, s) {
  ns <- length(s)
  hs <- numeric(ns)
  hs[1] <- 0.5 * log(s[2] / s[1])
  hs[ns] <- 0.5 * log(s[ns] / s[ns - 1])
  hs[2:(ns - 1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns - 2)]))
  
  n <- length(w)
  
  ws <- outer(w, s, "*")
  ws2 <- ws^2
  Hsuper <- matrix(exp(H), nrow = 2 * n, ncol = ns, byrow = TRUE) * rep(hs, each = 2 * n)
  
  DK <- rbind(ws2 / (1 + ws2), ws / (1 + ws2)) * Hsuper
  
  return(DK)
}

# Example usage
H <- runif(10)
w <- seq(1, 10, length.out = 10)
s <- seq(1, 10, length.out = 10)

DK_result <- kernelD(H, w, s)
print(DK_result)
