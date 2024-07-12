

LLS <- function(t, tau, Gexp) {
 
  
  n <- length(Gexp)
  
  res <- meshgrid(tau,t)
  S <- res$X
  T <- res$Y
  
  
  K <- (exp(-T/S))
  #rm(S, T)
  
  # gets (Gt/GtE - 1)^2, instead of  (Gt -  GtE)^2
  
  diag_matrix <- diag(1 / Gexp)

# Moltiplicare la matrice diagonale per K
  Kp <- diag_matrix %*% K
  
  Kp <- as.matrix(Kp)
  condKp <- pracma::cond(Kp)
  Gexp <- as.matrix(Gexp)
  g <- qr.solve(Kp, matrix(1, nrow = nrow(Gexp), ncol = 1))
  
  #rm(Kp)
  
  GtM <- K %*% g
  error <- sum((GtM / Gexp - 1)^2)
  
  # Risultati
  list(condKp = condKp, g = g, error = error)
  
}


#maxwell_modes <- MaxwellModes(z, t, Gt, par$prune)

MaxwellModes <- function(z, t, Gt, prune = 0) {

  N <- length(z)
  
  tau <- exp(z)
  n <- length(t)
  Gexp <- Gt
  
  prune <- ifelse(missing(prune), 0, prune)
  
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
    gneg <- g1[ineg]
    
    tau2 <- tau1
    tau2 <- tau2[-ineg]
    
    result2 <- LLS(t, tau2, Gexp)
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
  
  return(list(g = g, tau = tau, error = error, condKp = condKp))
  

}

