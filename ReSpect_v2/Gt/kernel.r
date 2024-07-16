# Function: kernel
#
# outputs the n*1 dimensional vector K(H)(t) which is comparable to Gexp = Gt
# modifying kernel for unevenly spaced s_i
#
# Input: H = substituted CRS,
#        t = n*1 vector contains times,
#        s = relaxation modes
#

kernel <- function(H, t, s) {
  
  ns <- length(s)
  hs <- numeric(ns)
  
  # Integration uses trapezoidal rule
  hs[1] <- 0.5 * log(s[2] / s[1])
  hs[ns] <- 0.5 * log(s[ns] / s[ns - 1])
  hs[2:(ns - 1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns - 2)]))
  
  res <- meshgrid(s,t)

  S <- res$X
  T <- res$Y 
  
  Kern = exp(-T/S)
  
  K <- Kern %*% (hs * exp(H))
  
  return(K)
}
