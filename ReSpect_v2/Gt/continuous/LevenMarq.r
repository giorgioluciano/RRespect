%# Function LevenMarq
#
# Purpose: Given a lambda, this function finds the H_lambda(s) that minimizes V(lambda)
#
#          V(lambda) := 1/n * ||Gst - kernel(H)||^2 +  lambda/nl * ||L H||^2
#
# Input  : lambda = regularization parameter,
#          Gst    = experimental data,
#          H      = guessed H,
#          t      = n*1 vector contains "times",
#          s      = relaxation modes,
#
# Output : H_lambda
#
#          Uses the Levenberg-Marquardt method to solve the nonlinear
#          minimization problem. Following the algorithm in the source:
#
# "A Brief Description of the Levenberg-Marquardt Algorithm Implemented by levmar
#  Manolis I. A. Lourakis"
#

LevenMarq <- function(lambda, Gst, H, t, s) {
  
  n <- length(t)
  ns <- length(s)
  nl <- ns - 2
  hs <- s[2] / s[1]
  
  tau <- 1e-3          # scaling parameter 1
  nu <- 2              # scaling parameter 2
  
  r <- rep(0, n + nl)
  Jr <- matrix(0, n + nl, ns)
  i <- 0
  pf <- sqrt(lambda)
  
  # L is a nl*ns tridiagonal matrix with 1, -2, and 1 on its diagonal.
  L <- diag(1, ns-1, ns) + diag(1, ns-1, ns) + diag(-2, ns, ns)
  L <- L[2:(nl+1), ]
  
  res <- GetResidualJacobian(pf, L, Gst, H, t, s)
  r <- res$r
  Jr <- res$Jr
  
  mu <- tau * max(diag(t(Jr) %*% Jr))  # scaling parameter 3
  
  # Doing the LM stuff
  ContinueCriteria <- TRUE
  
  while (ContinueCriteria && i < 5000) {
    
    i <- i + 1  # iteration count
    
    Delta <- solve(t(Jr) %*% Jr + mu * diag(ns), -t(Jr) %*% r)
    Hnew <- H + Delta
    
    res_new <- GetResidualJacobian(pf, L, Gst, Hnew, t, s)
    rnew <- res_new$r
    
    rho <- (norm(r, type="2")^2 - norm(rnew, type="2")^2) / (t(Delta) %*% (mu * Delta - t(Jr) %*% r))
    
    if (rho > 0) {
      
      H <- Hnew
      res <- GetResidualJacobian(pf, L, Gst, H, t, s)
      r <- res$r
      Jr <- res$Jr
      
      mu <- mu * max(1/3, 1 - (2 * rho - 1)^3)
      nu <- 2
      
      if ((norm(Jr %*% r, type="I") < 1e-6) || 
          ((norm(Delta, type="2") < 1e-6) || 
           (norm(Delta, type="2") < 1e-6 * norm(H, type="2")))) {
        
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

# HELPER FUNCTION: Gets Residuals r and Jacobian J
# If called with only one output variable, only returns residual

GetResidualJacobian <- function(pf, L, Gst, H, t, s) {
  
  n <- length(t)
  ns <- length(s)
  nl <- ns - 2
  
  r <- rep(0, n + nl)
  Jr <- matrix(0, n + nl, ns)
  
  # Get the residual vector first
  r[1:n] <- (1 - kernel(H, t, s) / Gst) / sqrt(n)
  r[(n+1):(n+nl)] <- pf * diff(H, differences = 2) / sqrt(nl)
  
  if (exists("Jr")) {
    
    Kmatrix <- (1 / Gst) * matrix(1, n, ns) / sqrt(n)
    Jr[1:n, ] <- -kernelD(H, t, s) * Kmatrix
    Jr[(n+1):(n+nl), ] <- pf * L / sqrt(nl)
    
  }
  
  return(list(r = r, Jr = Jr))
  
}

