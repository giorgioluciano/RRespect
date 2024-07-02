# Define the LevenMarq function as used previously

InitializeH <- function(Gexp, w, s) {
  # Initial guess
  H <- -5.0 + sin(pi * s)
  
  # Initial large lambda
  lambda <- 1e0
  Hlam <- LevenMarq(lambda, Gexp, H, w, s)
  
  # Successively improve the initial guess for low lambda
  lambda <- 1e-8
  H <- LevenMarq(lambda, Gexp, Hlam, w, s)
  
  return(H)
}

# Example usage
Gexp <- runif(20)  # Example experimental data
w <- seq(1, 10, length.out = 10)  # Example frequencies
s <- seq(1, 10, length.out = 10)  # Example relaxation modes

H <- InitializeH(Gexp, w, s)
print(H)
