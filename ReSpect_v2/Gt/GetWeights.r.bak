GetWeights <- function(H, t, s) {
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
 

  u=t  
  S <- outer(s, rep(1, length(t)))
  T <- outer(rep(1, length(s)), u)
  
  kern <- t(exp(-T/S))
  
  wij <-  kern %*% diag(hs * exp(H))
  K <- kern %*% (hs * exp(H))

  for (i in 1:(n)) {
    wij[i, ] <- wij[i, ] / K[i]
  }

  for (j in 1:ns) {
    wt[j] <- sum(wij[, j])
  }

  return(wt)
}

