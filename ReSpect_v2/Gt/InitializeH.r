# Define the LevenMarq function as used previously

InitializeH <- function(Gexp, t, s) {
  # Initial guess
  H <- -5.0 + sin(pi * s)
  
  # Initial large lambda
  lambda <- 1e0
  Hlam <- LevenMarq(lambda, Gexp, H, t, s)
  
  # Successively improve the initial guess for low lambda
  lambda <- 1e-8
  H <- LevenMarq(lambda, Gexp, Hlam, t, s)
  
  return(H)
}

