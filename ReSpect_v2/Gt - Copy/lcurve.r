#  
#  Function: lcurve(input)
# 
#  Input: Gexp = 2n*1 vector [Gp; Gpp],
#         Hgs = guessed H,
#         w = n*1 vector contains frequencies,
#         s = relaxation modes,
#         SmoothFac = Indirect way of controlling lambda_C, Set between -1 
#                    (lowest lambda explored) and 1 (highest lambda explored);
#                    When set to 0, using lambda_C determined from the
#                    L-curve,
# 
#  Output: lamC and 3 vectors of size npoints*1 contains a range of lambda, rho
#  and eta. "Elbow"  = lamC is estimated using a heuristic.
# 
#  Plot the L-curve. Can call the routine "corner.m" to find the elbow again

lcurve <- function(Gexp, Hgs, t, s, SmoothFac) {
  
  npoints <- 40
  lam_min <- 1e-12
  lam_max <- 1e+3
  
  hlam <- (lam_max/lam_min)^(1/(npoints-1))
  lam <- lam_min * hlam^(0:(npoints-1))
  
  eta <- numeric(npoints)
  rho <- numeric(npoints)
  
  # This step can be "parfor"ed if you have multicore
  
  for (i in 1:length(lam)) {
    lambda <- lam[i]
    H <- LevenMarq(lambda, Gexp, Hgs, t, s)
    rho[i] <- norm((1 - kernel(H, t, s) / Gexp), type = "2")
    eta[i] <- norm(diff(H, differences = 2), type = "2")
  }
  
  er <- rho/max(rho) + eta/max(eta)
  ermin <- min(er)
  i <- which.min(er)
  
  lamC <- lam[i]
  
  if (SmoothFac > 0) {
    lamC <- exp(log(lamC) + SmoothFac * (log(lam_max) - log(lamC)))
  } else if (SmoothFac < 0) {
    lamC <- exp(log(lamC) + SmoothFac * (log(lamC) - log(lam_min)))
  }
  
  return(list(lamC = lamC, lam = lam, rho = rho, eta = eta))
}
