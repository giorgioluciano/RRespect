

LLS <- function(w, tau, Gexp) {
  n <- length(Gexp) / 2
  X <- outer(tau, w, "*")
  ws <- X
  ws2 <- ws^2
  
  head_K <- (ws2 / (1 + ws2))
  tail_K <- (ws / (1 + ws2))
  
  K <- cbind(head_K, tail_K)
  K <-t(K)
  #check
  
  # Verify dimensions before multiplication
  if (nrow(K) != length(Gexp)) {
    stop("Dimension mismatch between K and Gexp")
  }
  
  
  Kp <- diag(1 / Gexp) %*% K
  
  condKp  = pracma::cond(Kp)

  ones_matrix <- matrix(1, nrow = length(Gexp))
  
  
  fit <- lm(ones_matrix ~ Kp - 1)  # -1 per evitare l'intercetta
 
  g=fit$coefficients
  
  
  GpM_nog <- t(ws2/(1+ws2))
  GpM <- GpM_nog  %*% g
  
  GppM_nog <- t(ws/(1+ws2))
  GppM <- GppM_nog  %*% g

  error <- sum((GpM / Gexp[1:n] - 1)^2 + (GppM / Gexp[(n + 1):(2 * n)] - 1)^2)
  
  
  list(g = g, error = error, condKp = condKp)
}


#maxwell_modes <- MaxwellModes(z, w, Gp, Gpp, par$prune)

MaxwellModes <- function(z, w, Gp, Gpp, prune = 0) {
  N <- length(z)
  tau <- exp(z)
  n <- length(w)
  Gexp <- c(Gp, Gpp)
  
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

