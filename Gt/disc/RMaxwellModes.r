# Function: MaxwellModes
# Solves the linear least squares problem to obtain the DRS
# Input:
#   z     = points distributed according to the density
#   t     = n*1 vector contains times
#   Gt    = n*1 vector contains G(t)
#   prune = Avoid modes with negative weights (1), or don't prune (0)
# Output:
#   g     = spectrum
#   tau   = array of relaxation times
#   error = relative error between the input data and the inferred G(t) from DRS
#   condKp = condition number
MaxwellModes <- function(z, t, Gt, prune = 0) {
  
  N <- length(z)
  tau <- exp(z)
  n <- length(t)
  Gexp <- Gt
  
  # Prune small negative weights g(i)
  if (!exists("prune")) {
    prune <- 0
  }
  
  if (!prune) {
    result <- LLS(t, tau, Gexp)
    g <- result$g
    error <- result$error
    condKp <- result$condKp
  } else {
    tau1 <- tau
    result1 <- LLS(t, tau1, Gexp)
    g1 <- result1$g
    error1 <- result1$error
    condKp1 <- result1$condKp
    
    # Run through negative "g"s and mark for deletion if the weight is small
    # Then repeat above calculation
    ineg <- which(g1 < 0)
    gneg <- g1[ineg]
    
    tau2 <- tau1
    tau2[ineg] <- NULL
    
    result2 <- LLS(t, tau2, Gexp)
    g2 <- result2$g
    error2 <- result2$error
    condKp2 <- result2$condKp
    
    if (condKp2 < condKp1 && (error2 / error1 - 1) < 0.05) {
      g <- g2
      error <- error2
      condKp <- condKp2
      tau <- tau2
    } else {
      g <- g1
      error <- error1
      condKp <- condKp1
      tau <- tau1
    }
  }
  
  return(list(g = g, tau = tau, error = error, condKp = condKp))
}

# Subfunction which does the actual LLS problem
LLS <- function(t, tau, Gexp) {
  n <- length(Gexp)
  K <- outer(tau, t, function(tau, t) exp(-t / tau))
  
  Kp <- solve(diag(1 / Gexp) %*% K)  # Solve the LLS problem
  
  condKp <- cond(Kp)  # Condition number of Kp
  g <- Kp %*% rep(1, n)  # Solve for g vector
  
  GtM <- K %*% g  # Compute G(t) from inferred spectrum
  error <- sum((GtM / Gexp - 1)^2)  # Compute relative error
  
  return(list(g = g, error = error, condKp = condKp))
}

# Example usage:
# Replace 'z', 't', and 'Gt' with appropriate vectors for your data
z <- log(c(0.1, 1, 10))  # Example points distributed according to density (logarithmically spaced)
t <- seq(0, 10, length.out = 100)  # Example vector of times
Gt <- sin(t) + rnorm(100, 0, 0.1)  # Example vector of G(t)

# Call the MaxwellModes function
result <- MaxwellModes(z, t, Gt, prune = 1)  # Adjust prune as needed (0 or 1)

# Print the results
print(result$g)      # Spectrum g
print(result$tau)    # Array of relaxation times tau
print(result$error)  # Relative error
print(result$condKp) # Condition number
