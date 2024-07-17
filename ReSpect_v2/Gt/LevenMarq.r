# Function LevenMarq
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
  
  K <- Kern * (hs %*% exp(H))
  
  return(K)
}


kernelD <- function(H, t, s) {
  
  ns <- length(s)
  hs <- numeric(ns)
  hs[1] <- 0.5 * log(s[2] / s[1])
  hs[ns] <- 0.5 * log(s[ns] / s[ns-1])
  hs[2:(ns-1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns-2)]))
  
  n <- length(t)
  
  res <- meshgrid(s,t)
  S =res$X
  Y =res$Y
  
  kern= exp(-T/S)
  
  
  Hsuper <- matrix(rep(hs * exp(H), each = n), nrow = n, ncol = ns, byrow = TRUE)
  
  DK <- kern * Hsuper
  
  return(DK)
}


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
  
  L <- diag(-2, ns)
  L[row(L) == col(L) - 1] <- 1
  L[row(L) == col(L) + 1] <- 1
  L_sub <- L[2:(nl+1), ]
  L <- L_sub
  
  
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
    
    if (!is.na(rho) && rho > 0) {
      
      H <- Hnew
      res <- GetResidualJacobian(pf, L, Gst, H, t, s)
      r <- res$r
      Jr <- res$Jr
      
      mu <- mu * max(1/3, 1 - (2 * rho - 1)^3)
      nu <- 2
      
      norm_Jr_r_inf <- norm(t(Jr) %*% r, type = "I")
      norm_Delta_2 <- norm(Delta, type = "2")
      norm_H_2 <- norm(H, type = "2")
      
      if ( (norm_Jr_r_inf < 1e-6) || 
           ((norm_Delta_2 < 1e-6) || 
            (norm_Delta_2 < 1e-6 * norm_H_2)) ) {
        
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


