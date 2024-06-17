# Function: lcurve
# Output: lamC, lam, rho, eta
# Input: Gexp = 2n*1 vector [Gp; Gpp],
#        Hgs = guessed H,
#        w = n*1 vector contains frequencies,
#        s = relaxation modes,
#        SmoothFac = Indirect way of controlling lambda_C, Set between -1 
#                   (lowest lambda explored) and 1 (highest lambda explored);
#                   When set to 0, using lambda_C determined from the
#                   L-curve,
lcurve <- function(Gexp, Hgs, w, s, SmoothFac) {
  npoints <- 40
  lam_min <- 1e-12
  lam_max <- 1e+3
  
  hlam <- (lam_max / lam_min)^(1 / (npoints - 1))
  lam <- lam_min * hlam^(0:(npoints - 1))
  rho <- numeric(npoints)
  eta <- numeric(npoints)
  
  # Compute rho and eta for each lambda
  for (i in seq_along(lam)) {
    lambda <- lam[i]
    H <- LevenMarq(lambda, Gexp, Hgs, w, s)
    K <- kernel(H, w, s)
    rho[i] <- norm(1 - K / Gexp)
    eta[i] <- norm(diff(H, 2))
  }
  
  # Find lamC based on the L-curve heuristic
  er <- rho / max(rho) + eta / max(eta)
  i <- which.min(er)
  lamC <- lam[i]
  
  # Adjust lamC based on SmoothFac
  if (SmoothFac > 0) {
    lamC <- exp(log(lamC) + SmoothFac * (log(lam_max) - log(lamC)))
  } else if (SmoothFac < 0) {
    lamC <- exp(log(lamC) + SmoothFac * (log(lamC) - log(lam_min)))
  }
  
  return(list(lamC = lamC, lam = lam, rho = rho, eta = eta))
}

# Example usage:
# Ensure you have LevenMarq and kernel functions defined properly in your R environment
# Define Gexp, Hgs, w, s, and SmoothFac appropriately before calling lcurve

# par <- SetParameters()  # If using SetParameters function
# Gexp <- GetExpData(par$GstFile)
# Hgs <- InitializeH(Gexp, w, s)
# lcurve_result <- lcurve(Gexp, Hgs, w, s, par$SmFacLam)
# lamC <- lcurve_result$lamC
# lam <- lcurve_result$lam
# rho <- lcurve_result$rho
# eta <- lcurve_result$eta

