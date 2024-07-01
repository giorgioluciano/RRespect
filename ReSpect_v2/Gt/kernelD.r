# Function: kernelD
#
# outputs the n*ns dimensional vector DK(H)(t)
# approximates dK/dHj
#
# Input: H = substituted CRS,
#        t = n*1 vector of times,
#        s = relaxation modes,
#
# Output: DK = Jacobian of H
#

kernelD <- function(H, t, s) {
  
  ns <- length(s)
  hs <- numeric(ns)
  hs[1] <- 0.5 * log(s[2] / s[1])
  hs[ns] <- 0.5 * log(s[ns] / s[ns-1])
  hs[2:(ns-1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns-2)]))
  
  n <- length(t)
  
  kern <- outer(t, s, function(T, S) exp(-T / S))
  
  Hsuper <- matrix(rep(hs * exp(H), each = n), nrow = n, ncol = ns, byrow = TRUE)
  
  DK <- kern * Hsuper
  
  return(DK)
}
