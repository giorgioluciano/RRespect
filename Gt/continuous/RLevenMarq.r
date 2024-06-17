# Function: LevenMarq
# Purpose: Given a lambda, this function finds the H_lambda(s) that minimizes V(lambda)
#          V(lambda) := 1/n * ||Gst - kernel(H)||^2 +  lambda/nl * ||L H||^2
# Input  : lambda = regularization parameter,
#          Gst    = experimental data,
#          H      = guessed H,
#          t      = n*1 vector contains "times",
#          s      = relaxation modes,
# Output : Hlambda
# Uses the Levenberg-Marquardt method to solve the nonlinear
# minimization problem.
LevenMarq <- function(lambda, Gst, H, t, s) {
  
  n <- length(t)
  ns <- length(s)
  nl <- ns - 2
  hs <- s[2] / s[1]
  
  tau <- 1e-3    # scaling parameter 1
  nu <- 2        # scaling parameter 2
  
  r <- rep(0, n + nl)
  Jr <- matrix(0, n + nl, ns)
  i <- 0
  pf <- sqrt(lambda)
  
  # L is a nl*ns tridiagonal matrix with 1, -2, and 1 on its diagonal.
  L <- diag(1, nl) + diag(1, nl - 1, k = 1) + diag(-2, ns)
  L <- L[-1, ]  # Remove the first row
  
  # Initial calculation of residual r and Jacobian Jr
  result <- GetResidualJacobian(pf, L, Gst, H, t, s)
  r <- result$r
  Jr <- result$Jr
  
  mu <- tau * max(diag(t(Jr) %*% Jr))  # scaling parameter 3
  
  ContinueCriteria <- TRUE
  
  while (ContinueCriteria && i < 5000) {
    i <- i + 1  # iteration count
    
    Delta <- solve(t(Jr) %*% Jr + mu * diag(ns)) %*% (-t(Jr) %*% r)
    Hnew <- H + Delta
    
    result <- GetResidualJacobian(pf, L, Gst, Hnew, t, s)
    rnew <- result$r
    
    rho <- (sum(r^2) - sum(rnew^2)) / (t(Delta) %*% (mu * Delta - t(Jr) %*% r))
    
    if (rho > 0) {
      H <- Hnew
      result <- GetResidualJacobian(pf, L, Gst, H, t, s)
      r <- result$r
      Jr <- result$Jr
      
      mu <- mu * max(1/3, 1 - (2 * rho - 1)^3)
      nu <- 2
      
      # Check stopping criteria
      if (max(abs(t(Jr) %*% r)) < 1e-6 || 
          max(abs(Delta)) < 1e-6 || 
          max(abs(Delta)) < 1e-6 * max(abs(H))) {
        ContinueCriteria <- FALSE
      }
      
    } else {
      mu <- mu * nu
      nu <- 2 * nu
    }
  }
  
  Hlambda <- H
  
  return(Hlambda)
}

# Helper Function: GetResidualJacobian
# Gets Residuals r and Jacobian J
# If called with only one output variable, only returns residual
GetResidualJacobian <- function(pf, L, Gst, H, t, s) {
  
  n <- length(t)
  ns <- length(s)
  nl <- ns - 2
  
  r <- rep(0, n + nl)
  Jr <- matrix(0, n + nl, ns)
  
  # Get the residual vector r
  r[1:n] <- (1 - kernel(H, t, s) / Gst) / sqrt(n)
  r[(n + 1):(n + nl)] <- pf * diff(H, differences = 2) / sqrt(nl)
  
  # Furnish the Jacobian Jr
  if (!missing(pf)) {
    Kmatrix <- matrix(1 / Gst, n, ns) / sqrt(n)
    Jr[1:n, ] <- -kernelD(H, t, s) * Kmatrix
    Jr[(n + 1):(n + nl), ] <- pf * L / sqrt(nl)
  }
  
  return(list(r = r, Jr = Jr))
}

# Example usage:
# Replace 'lambda', 'Gst', 'H', 't', and 's' with appropriate vectors/values for your data
lambda <- 0.1    # Example regularization parameter
Gst <- runif(100)  # Example experimental data vector
H <- runif(5)     # Example guessed H vector
t <- seq(0, 10, length.out = 100)  # Example time vector
s <- seq(1, 5, by = 1)             # Example relaxation modes

# Call the LevenMarq function
Hlambda <- LevenMarq(lambda, Gst, H, t, s)

