
# Function LevenMarq(input)
#
# Purpose: Given a lambda, this function finds the H_lambda(s) that minimizes V(lambda)

#          V(lambda) := 1/n * ||Gst - kernel(H)||^2 + lambda/nl * ||L H||^2

# Input  : lambda = regularization parameter ,
#          Gst = experimental data,
#          H = guessed H,
#          w = n*1 vector contains frequencies,
#          s = relaxation modes,

# Output : H_lambda
 
#          Uses the Levenberg-Marquardt method to solve the nonlinear
#          minimization problem. Following the algorithm in the source:

# "A Brief Description of the Levenberg-Marquardt Algorithm Implemened by levmar
#  Manolis I. A. Lourakis"


library(Matrix)

GetResidualJacobian <- function(pf, L, Gst, H, t, s) {
  n <- length(t)
  ns <- length(s)
  nl <- ns - 2
  
  r <- numeric(n + nl)
  Jr <- matrix(0, nrow = n + nl, ncol = ns)
  
  # Get the residual vector first
  # r = vector of size (2n+nl,1)
  
  head_r <-(1 - kernel(H, t, s) / Gst) / sqrt(n)
  tail_r <- pf * diff(diff(H)) / sqrt(nl)  # second derivative
  
  res <- c(head_r,tail_r)
  
  # Furnish the Jacobian Jr
  # (2n+nl)*ns matrix
  
  
  Kmatrix <- matrix(1/Gst, nrow = n, ncol = ns) / sqrt(n)
  
  head_Jr <- -kernelD(H, t, s) * Kmatrix
  tail_Jr <- pf * L / sqrt(nl)
  
  Jr <- rbind(head_Jr,tail_Jr)
  return(list(r = res, Jr = Jr))
}


kernel <- function(H, t, s) {
  ns <- length(t)
  hs <- numeric(ns)
  
  # Uses trapezoidal rule for integration
  hs[1] <- 0.5 * log(s[2] / s[1])
  hs[ns] <- 0.5 * log(s[ns] / s[ns - 1])
  hs[2:(ns - 1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns - 2)]))
  
  # Create meshgrid equivalent
  
  res <- meshgrid(s,t)
  S <- res$X
  T <- res$Y 
  
  kern <- t(exp(-T/S))

  #rm(S, T)
  
  K <- kern %*% (hs * exp(H))
  
  return(K)
}


kernelD <- function(H, t, s) {
  ns <- length(s)
  hs <- numeric(ns)
  hs[1] <- 0.5 * log(s[2] / s[1])
  hs[ns] <- 0.5 * log(s[ns] / s[ns - 1])
  hs[2:(ns - 1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns - 2)]))
  
  n <- length(t)
  
  res <- meshgrid(s,t)
  S <- res$X
  T <- res$Y 
  
  
  kern <- t(exp(-T/S))
  
  Hsuper <- matrix(exp(H), nrow = n, ncol = ns, byrow = TRUE) * rep(hs, each =  n)
  
  DK <- kern * Hsuper
  
  return(DK)
}



LevenMarq <- function(lambda, Gst, H, t, s) {
  n <- length(t)
  ns <- length(s)
  nl <- ns - 2
  hs <- s[2] / s[1]
  
  tau <- 1e-3
  nu <- 2
  
  r <- numeric(n + nl)
  Jr <- matrix(0, n + nl, ns)
  i <- 0
  pf <- sqrt(lambda)
  
  
  # Create the tridiagonal matrix L
  
  L <- diag(-2, ns)
  L[row(L) == col(L) - 1] <- 1
  L[row(L) == col(L) + 1] <- 1
  L_sub <- L[2:(nl+1), ]
  L <- L_sub
  
  residuals_and_jacobian <- GetResidualJacobian(pf, L, Gst, H, t, s)
  
  r <- residuals_and_jacobian$r
  Jr <- residuals_and_jacobian$Jr
  
  
  mu <- tau * max(diag(crossprod(Jr)))
  
  ContinueCriteria <- TRUE
  
  while (ContinueCriteria && i < 5000) {
    i <- i + 1
    
    Delta <- qr.solve(crossprod(Jr) + mu * diag(ns), -crossprod(Jr, r))
    Hnew <- H + Delta
    
    new_residuals <- GetResidualJacobian(pf, L, Gst, Hnew, t, s)$r
    rho <- (sum(r^2) - sum(new_residuals^2)) / (crossprod(Delta, mu * Delta - crossprod(Jr, r)))
    
    if (rho > 0) {
      H <- Hnew
      residuals_and_jacobian <- GetResidualJacobian(pf, L, Gst, H, t, s)
     
	 r <- residuals_and_jacobian$r
      Jr <- residuals_and_jacobian$Jr
      
      mu <- mu * max(1 / 3, 1 - (2 * rho - 1)^3)
      nu <- 2
      
      if ((norm(Delta, type = "2") < 1e-6 && norm(crossprod(Jr, r), type = "I") < 1e-6) ||
          (norm(Delta, type = "2") < 1e-6 * norm(H, type = "2"))) {
        ContinueCriteria <- FALSE
      }
    } else {
      mu <- mu * nu
      nu <- 2 * nu
    }
  }
  
  return(H)
}



