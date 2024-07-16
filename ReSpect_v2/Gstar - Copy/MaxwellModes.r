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
MaxwellModes <- function(z, w, Gp, Gpp, prune = FALSE) {
  
  N <- length(z)
  tau <- exp(z)
  n <- length(w)
  Gexp <-  c(Gp, Gpp)
  
  # Prune small negative weights g(i)
  if (!prune) {
    result <- LLS(w, tau, Gexp)
    g <- result$g
    error <- result$error
    condKp <- result$condKp
    
  } else {
    tau1 <- tau
    result1 <- LLS(w, tau1, Gexp)
    
    g1 <- result1$g
    error1 <- result1$error
    condKp1 <- result1$condKp
    
    # Identify negative "g"s and mark for deletion if weight is small
    ineg <- which(g1 < 0)
    gneg <- g1[ineg]
    
    tau2 <- tau1
    
    tau2 <- tau2[-ineg]
    
    result2 <- LLS(w, tau2, Gexp)
    g2 <- result2$g
    error2 <- result2$error
    condKp2 <- result2$condKp
    
    if (condKp2 < condKp1 && (error2/error1 - 1) < 0.05) {
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
LLS <- function(w, tau, Gexp) {
  
  n <- nrow(Gexp) / 2
  resi <- meshgrid(tau,w)
  X <- resi$X
  Y <- resi$Y
  
 
  ws <- X * Y
  ws2 <- ws^2
  

  K <- rbind(ws2 / (1 + ws2), ws / (1 + ws2))
  
  diag_matrix <- diag(1 / Gexp)
  
  Kp <- diag_matrix %*% K
  
  Kp <- as.matrix(Kp)
  condKp <- pracma::cond(Kp)
  Gexp <- as.matrix(Gexp)
  g <- qr.solve(Kp, matrix(1, nrow = nrow(Gexp), ncol = 1))
  
  condKp <- cond(Kp)
  
  ones_vector <- rep(1, length(Gexp))
  g <- qr.solve(Kp, ones_vector)
    
  GpM <- (ws2 / (1 + ws2)) %*% g
  GppM <- (ws / (1 + ws2)) %*% g
  
  GpM_err <- (GpM / Gexp[1:length(GpM), ] - 1)^2
  GppM_err <- GppM / (Gexp[(length(GppM) + 1):(2 * length(GppM)), ] - 1)^2
  
  error <- sum( GpM_err + GppM_err)
  
  return(list(g = g, error = error, condKp = condKp)) 

}
