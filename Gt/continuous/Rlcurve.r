# Function: lcurve
# Output: lamC, lam, rho, eta
# Input: Gexp = n*1 vector [Gt],
#        Hgs  = guessed H,
#        t    = n*1 vector contains times,
#        s    = relaxation modes,
#        SmoothFac = Indirect way of controlling lambda_C,
#                    Set between -1 (lowest lambda explored) and 1 (highest lambda explored);
#                    When set to 0, using lambda_C determined from the L-curve,
lcurve <- function(Gexp, Hgs, t, s, SmoothFac) {
  
  npoints <- 40
  lam_min <- 1e-12
  lam_max <- 1e+3
  
  hlam <- (lam_max / lam_min)^(1 / (npoints - 1))
  lam <- lam_min * hlam^(0:(npoints - 1))
  
  eta <- numeric(npoints)
  rho <- numeric(npoints)
  
  # Calculation of rho and eta
  for (i in 1:length(lam)) {
    lambda <- lam[i]
    H <- LevenMarq(lambda, Gexp, Hgs, t, s)  # Assuming LevenMarq is defined
    rho[i] <- norm((1 - kernel(H, t, s) / Gexp))
    eta[i] <- norm(diff(H, differences = 2))
  }
  
  # Finding lamC using a heuristic based on rho and eta
  er <- rho / max(rho) + eta / max(eta)
  i <- which.min(er)
  lamC <- lam[i]
  
  # Adjusting lamC based on SmoothFac
  if (SmoothFac > 0) {
    lamC <- exp(log(lamC) + SmoothFac * (log(lam_max) - log(lamC)))
  } else if (SmoothFac < 0) {
    lamC <- exp(log(lamC) + SmoothFac * (log(lamC) - log(lam_min)))
  }
  
  return(list(lamC = lamC, lam = lam, rho = rho, eta = eta))
}

# Example usage:
# Replace 'Gexp', 'Hgs', 't', 's', and 'SmoothFac' with appropriate vectors/values for your data
Gexp <- runif(100)  # Example Gexp vector
Hgs <- runif(5)     # Example guessed H vector
t <- seq(0, 10, length.out = 100)  # Example time vector
s <- seq(1, 5, by = 1)             # Example relaxation modes
SmoothFac <- 0.2    # Example smoothness factor

# Call the lcurve function
result <- lcurve(Gexp, Hgs, t, s, SmoothFac)

# Access the outputs
lamC <- result$lamC
lam <- result$lam
rho <- result$rho
eta <- result$eta



