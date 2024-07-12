#
# Function: MaxwellModes(input)
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

LLS <- function(w, tau, Gexp) {
  n <- length(Gexp) / 2
  
  res <-  meshgrid(tau,w)
  X <- res$X
  Y <- res$ Y
  ws <- X  *Y 
  ws2 <- ws^2
  
  head_K <- (ws2 / (1 + ws2))
  tail_K <- (ws / (1 + ws2))
 
   K <- cbind(head_K,tail_K)
  
  Kp <- diag(1 / Gexp) %*% K
  Kp <- as.matrix(Kp)
  condKp <- pracma::cond(Kp)
  
  Gexp <- as.matrix(Gexp)
  g <- qr.solve(Kp, matrix(1, nrow = nrow(Gexp), ncol = 1))
  
  K <- exp(-T / S)
  
  GpM <- head_K *g
  GppM <- tail_k*g
  
  error <- sum((GpM / Gexp[1:n] - 1)^2 + (GppM / Gexp[(n+1):(2*n)] - 1)^2)
  
  list(g = g, error = error, condKp = condKp)
}

#maxwell_modes <- MaxwellModes(z, w, Gp, Gpp, par$prune)

MaxwellModes <- function(z, t, Gt,  prune = 0) {
  N <- length(z)
  tau <- exp(z)
  n <- length(w)
  Gexp <- Gt
  
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
    
    ineg <- which(g1 < 0)
    gneg <- g1[ineg]
    
    tau2 <- tau1
    tau2 <- tau2[-ineg]
    
    result2 <- LLS(w, tau2, Gexp)
  
    g2 <- result2$g
    error2 <- result2$error
    condKp2 <- result2$condKp
    
    if (condKp2 < condKp1 && (error2 / error1 - 1) < 0.05) {
      error <- error2
      condKp <- condKp2
      g <- g2
      tau <- tau2
    } else {
      error <- error1
      condKp <- condKp1
      g <- g1
      tau <- tau1
    }
  }
  
  list(g = g, tau = tau, error = error, condKp = condKp)
}

