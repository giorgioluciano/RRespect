# Function: LevenMarq
# Purpose: Given a lambda, this function finds the H_lambda(s) that minimizes V(lambda)
#          V(lambda) := 1/n * ||Gst - kernel(H)||^2 + lambda/nl * ||L H||^2
# Input: lambda = regularization parameter ,
#        Gst = experimental data,
#        H = guessed H,
#        w = n*1 vector contains frequencies,
#        s = relaxation modes,
# Output: Hlambda
#         Uses the Levenberg-Marquardt method to solve the nonlinear
#         minimization problem.
LevenMarq <- function(lambda, Gst, H, w, s) {
  n <- length(w)
  ns <- length(s)
  nl <- ns - 2
  hs <- s[2] / s[1]
  
  tau <- 1e-3          # scaling parameter 1
  nu <- 2              # scaling parameter 2
  
  r <- numeric(2 * n + nl)
  Jr <- matrix(0, nrow = 2 * n + nl, ncol = ns)
  i <- 0
  pf <- sqrt(lambda)
  
  # L is a nl*ns tridiagonal matrix with 1, -2, and 1 on its diagonal
  L <- diag(1, ns - 1) + diag(1, ns - 1, 1) + diag(-2, ns)
  L <- L[2:(nl + 1), ]
  
  # Helper function to get residuals and Jacobian
  GetResidualJacobian <- function(pf, L, Gst, H, w, s) {
    r <- numeric(2 * n + nl)
    Jr <- matrix(0, nrow = 2 * n + nl, ncol = ns)
    
    # Get the residual vector r
    r[1:2 * n] <- (1 - kernel(H, w, s) / Gst) / sqrt(n)  # the Gp and Gpp residuals
    r[(2 * n + 1):(2 * n + nl)] <- pf * diff(H, 2) / sqrt(nl)  # second derivative residuals
    
    # Furnish the Jacobian Jr
    if (nrow(Jr) > 0) {
      Kmatrix <- (1 / Gst) * matrix(1, nrow = 1, ncol = ns) / sqrt(n)
      Jr[1:(2 * n), 1:ns] <- -kernelD(H, w, s) * Kmatrix
      Jr[(2 * n + 1):(2 * n + nl), 1:ns] <- pf * L / sqrt(nl)
    }
    
    return(list(r = r, Jr = Jr))
  }
  
  # Initial residual and Jacobian
  res <- GetResidualJacobian(pf, L, Gst, H, w, s)
  r <- res$r
  Jr <- res$Jr
  
  # Initial scaling parameter for LM
  mu <- tau * max(diag(t(Jr) %*% Jr))  # scaling parameter 3
  
  # Levenberg-Marquardt iteration
  while (i < 5000) {
    i <- i + 1  # iteration count
    
    # Compute the update Delta
    Delta <- solve(t(Jr) %*% Jr + mu * diag(ns)) %*% (-t(Jr) %*% r)
    Hnew <- H + Delta
    
    # Compute new residuals and Jacobian
    res <- GetResidualJacobian(pf, L, Gst, Hnew, w, s)
    rnew <- res$r
    
    # Compute rho
    rho <- (norm(r)^2 - norm(rnew)^2) / (t(Delta) %*% (mu * Delta - t(Jr) %*% r))
    
    if (rho > 0) {
      H <- Hnew
      r <- rnew
      Jr <- res$Jr
      
      mu <- mu * max(1/3, 1 - (2 * rho - 1)^3)
      nu <- 2
      
      # Stopping criteria
      if (norm(Delta) < 1e-6 && norm(t(Jr) %*% r, "inf") < 1e-6) {
        break
      }
    } else {
      mu <- mu * nu
      nu <- 2 * nu
    }
  }
  
  Hlambda <- H
  return(Hlambda)
}

# Example usage:
# Ensure you have GetResidualJacobian, kernel, and kernelD functions defined properly in your R environment
# Define lambda, Gst, H, w, and s appropriately before calling LevenMarq

# Hlambda <- LevenMarq(lambda, Gst, H, w, s)

