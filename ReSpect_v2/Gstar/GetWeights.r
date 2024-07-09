GetWeights <- function(H, w, s) {
  # Function: GetWeights(input)
  #
  # Finds the weight of "each" mode by taking a weighted average of its contribution
  # to Gp and Gpp
  #
  # Input: H = substituted CRS
  #        w = n*1 vector contains frequencies
  #        s = relaxation modes
  #
  # Output: wt = weight of each mode
  #

  ns <- length(s)
  n <- length(t)

  hs <- numeric(ns)
  wt <- numeric(ns)
  
  hs[1] <- 0.5 * log(s[2] / s[1])
  hs[ns] <- 0.5 * log(s[ns] / s[ns-1])
  hs[2:(ns-1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns-2)]))

 
  S <- outer(s, t, FUN = function(x, y) x)
  T <- outer(s, t, FUN = function(x, y) y)
  
  # Calcolo del kernel
  kern <- exp(-T / S)
  
  # Pulizia delle variabili S e T (non necessaria in R)
  rm(S, T)
  
  # Calcolo di wij e K
  wij <- kern %*% diag(hs * exp(H))
  K <- kern %*% (hs * exp(H))
  
  
  for (i in 1:(2 * n)) {
    wij[i, ] <- wij[i, ] / K[i]
  }

  for (j in 1:ns) {
    wt[j] <- sum(wij[, j])
  }

  return(wt)
}

