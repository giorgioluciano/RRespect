
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

LevenMarq <- function(lambda, Gst, H, w, s) {
  n <- length(w)
  ns <- length(s)
  nl <- ns - 2
  hs <- s[2] / s[1]
  
  tau <- 1e-3
  nu <- 2
  r <- numeric(2 * n + nl)
  Jr <- matrix(0, 2 * n + nl, ns)
  i <- 0
  pf <- sqrt(lambda)
  
  # Create the tridiagonal matrix L
  L <- bandSparse(nl, ns, c(-1, 0, 1), list(rep(1, nl), rep(-2, ns), rep(1, nl)))
  
  residuals_and_jacobian <- GetResidualJacobian(pf, L, Gst, H, w, s)
  r <- residuals_and_jacobian$residuals
  Jr <- residuals_and_jacobian$jacobian
  mu <- tau * max(diag(crossprod(Jr)))
  
  ContinueCriteria <- TRUE
  
  while (ContinueCriteria && i < 5000) {
    i <- i + 1
    
    Delta <- solve(crossprod(Jr) + mu * diag(ns), -crossprod(Jr, r))
    Hnew <- H + Delta
    
    new_residuals <- GetResidualJacobian(pf, L, Gst, Hnew, w, s)$residuals
    rho <- (sum(r^2) - sum(new_residuals^2)) / (crossprod(Delta, mu * Delta - crossprod(Jr, r)))
    
    if (rho > 0) {
      H <- Hnew
      residuals_and_jacobian <- GetResidualJacobian(pf, L, Gst, H, w, s)
      r <- residuals_and_jacobian$residuals
      Jr <- residuals_and_jacobian$jacobian
      
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


kernel <- function(H, w, s) {
  ns <- length(s)
  hs <- numeric(ns)
  
  # Uses trapezoidal rule for integration
  hs[1] <- 0.5 * log(s[2] / s[1])
  hs[ns] <- 0.5 * log(s[ns] / s[ns - 1])
  
  hs[2:(ns - 1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns - 2)]))
  
  # Create meshgrid equivalent
  ws <- outer(w, s, "*")
  ws2 <- ws^2
  
  # Calculate K
  K <- c((ws2 / (1 + ws2)) %*% (hs * exp(H)), (ws / (1 + ws2)) %*% (hs * exp(H)))
  
  return(K)
}


kernelD <- function(H, w, s) {
  ns <- length(s)
  hs <- numeric(ns)
  hs[1] <- 0.5 * log(s[2] / s[1])
  hs[ns] <- 0.5 * log(s[ns] / s[ns - 1])
  hs[2:(ns - 1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns - 2)]))
  
  n <- length(w)
  
  ws <- outer(w, s, "*")
  ws2 <- ws^2
  Hsuper <- matrix(exp(H), nrow = 2 * n, ncol = ns, byrow = TRUE) * rep(hs, each = 2 * n)
  
  DK <- rbind(ws2 / (1 + ws2), ws / (1 + ws2)) * Hsuper
  
  return(DK)
}

GetResidualJacobian <- function(pf, L, Gst, H, w, s) {
  n <- length(w)
  ns <- length(s)
  nl <- ns - 2
  r <- numeric(2*n + nl)
  Jr <- matrix(0, nrow = 2*n + nl, ncol = ns)
  
  # Get the residual vector first
  # r = vector of size (2n+nl,1)
  
  # Assumendo che kernel sia una funzione definita dall'utente
  r[1:(2*n)] <- (1 - kernel(H, w, s) / Gst) / sqrt(n)  # the Gp and Gpp
  
  # In R, diff() restituisce un vettore di lunghezza n-1, quindi usiamo diff() due volte
  r[(2*n+1):(2*n+nl)] <- pf * diff(diff(H)) / sqrt(nl)  # second derivative
  
  # Furnish the Jacobian Jr
  # (2n+nl)*ns matrix
  
  # In R, non abbiamo bisogno di controllare il numero di output
  Kmatrix <- matrix(1/Gst, nrow = 2*n, ncol = ns) / sqrt(n)
  
  # Assumendo che kernelD sia una funzione definita dall'utente
  Jr[1:(2*n), 1:ns] <- -kernelD(H, w, s) * Kmatrix
  Jr[(2*n+1):(2*n+nl), 1:ns] <- pf * L / sqrt(nl)
  
  # Restituisci una lista con r e Jr
  return(list(r = r, Jr = Jr))
}



