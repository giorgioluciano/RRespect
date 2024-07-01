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

GetResidualJacobian <- function(pf, L, Gst, H, w, s) {
  n <- length(w)
  ns <- length(s)
  nl <- ns - 2
  r <- numeric(2 * n + nl)
  Jr <- matrix(0, 2 * n + nl, ns)
  
  kernel_values <- kernel(H, w, s)
  r[1:2 * n] <- (1 - kernel_values / Gst) / sqrt(n)
  r[2 * n + 1:2 * n + nl] <- pf * diff(H, differences = 2) / sqrt(nl)
  
  if (nargout == 2) {
    Kmatrix <- (1 / Gst) %*% t(rep(1, ns)) / sqrt(n)
    Jr[1:2 * n, 1:ns] <- -kernelD(H, w, s) * Kmatrix
    Jr[2 * n + 1:2 * n + nl, 1:ns] <- pf * L / sqrt(nl)
  }
  
  return(list(residuals = r, jacobian = Jr))
}

kernel <- function(H, w, s) {
  # Dummy kernel function to be implemented
  return(rep(1, length(w))) # Replace with actual implementation
}

kernelD <- function(H, w, s) {
  # Dummy kernel derivative function to be implemented
  return(matrix(1, nrow = length(w), ncol = length(H))) # Replace with actual implementation
}

lcurve <- function(Gexp, Hgs, w, s, SmoothFac) {
  npoints <- 40
  
  lam_min <- 1e-12
  lam_max <- 1e+3
  
  hlam <- (lam_max / lam_min)^(1 / (npoints - 1))
  lam <- lam_min * hlam^seq(0, npoints - 1)
  
  eta <- numeric(npoints)
  rho <- numeric(npoints)
  
  for (i in 1:length(lam)) {
    lambda <- lam[i]
    H <- LevenMarq(lambda, Gexp, Hgs, w, s)
    rho[i] <- norm((1 - kernel(H, w, s) / Gexp))
    eta[i] <- norm(diff(H, differences = 2))
  }
  
  er <- rho / max(rho) + eta / max(eta)
  lamC <- lam[which.min(er)]
  
  if (SmoothFac > 0) {
    lamC <- exp(log(lamC) + SmoothFac * (log(lam_max) - log(lamC)))
  } else if (SmoothFac < 0) {
    lamC <- exp(log(lamC) + SmoothFac * (log(lamC) - log(lam_min)))
  }
  
  return(list(lamC = lamC, lam = lam, rho = rho, eta = eta))
}

# Example usage
Gexp <- runif(20)
Hgs <- runif(10)
w <- seq(1, 10, length.out = 10)
s <- seq(1, 10, length.out = 10)
SmoothFac <- 0

result <- lcurve(Gexp, Hgs, w, s, SmoothFac)
print(result)
