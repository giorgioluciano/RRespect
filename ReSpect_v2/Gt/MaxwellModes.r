# Function: MaxwellModes
#
# Solves the linear least squares problem to obtain the DRS
#
# Input: z = points distributed according to the density,
#        t  = n*1 vector contains times,
#        Gt = n*1 vector contains G(t),
#
#        Prune = Avoid modes with -ve weights (=1), or don't care (0) 
#
# Output: g, tau = spectrum  (array)
#         error = relative error between the input data and the G(t) inferred from the DRS
#         condKp = condition number
#

MaxwellModes <- function(z, t, Gt, prune = 0) {
  
  N <- length(z)
  tau <- exp(z)
  n <- length(t)
  Gexp <- Gt
  
  if (prune == 0) {
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
    
    ineg <- which(g1 < 0)
    tau2 <- tau1[-ineg]
    
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
  K <- outer(t, tau, function(T, S) exp(-T / S))
  
  Kp <- diag(1 / Gexp) %*% K
  condKp <- kappa(Kp)
  g <- solve(Kp, rep(1, n))
  
  GtM <- K %*% g
  error <- sum((GtM / Gexp - 1)^2)
  
  return(list(g = g, error = error, condKp = condKp))
}
