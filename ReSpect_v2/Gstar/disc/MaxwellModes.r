library(MASS)  # For ginv function (Moore-Penrose pseudoinverse)

MaxwellModes <- function(z, w, Gp, Gpp, prune = 0) {
  # Function: MaxwellModes(input)
  #
  # Solves the linear least squares problem to obtain the DRS
  #
  # Input: z = points distributed according to the density,
  #        w = n*1 vector contains frequencies,
  #        Gp, Gpp
  #        Prune = Avoid modes with -ve weights (=1), or don't care (0) 
  #
  # Output: g, tau = spectrum 
  #         error = relative error between the input data and the G*(w) inferred from the DRS
  #         condKp = condition number
  #

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

LLS <- function(w, tau, Gexp) {
  n <- length(Gexp) / 2
  X <- outer(tau, w, "*")
  ws <- X
  ws2 <- ws^2

  K <- rbind((ws2 / (1 + ws2)), (ws / (1 + ws2)))
  Kp <- diag(1 / Gexp) %*% K

  condKp <- kappa(Kp)
  g <- ginv(Kp) %*% rep(1, length(Gexp))

  GpM <- (ws2 / (1 + ws2)) %*% g
  GppM <- (ws / (1 + ws2)) %*% g
  error <- sum((GpM / Gexp[1:n] - 1)^2 + (GppM / Gexp[(n + 1):(2 * n)] - 1)^2)
  
  list(g = g, error = error, condKp = condKp)
}

# Example usage
# Define example input values
z <- c(1, 2, 3)
w <- seq(0.1, 10, length.out = 100)
Gp <- runif(100, min = 0, max = 10)
Gpp <- runif(100, min = 0, max = 10)

# Call the function
result <- MaxwellModes(z, w, Gp, Gpp, prune = 1)
print(result)
