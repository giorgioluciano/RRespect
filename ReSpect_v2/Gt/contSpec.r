# Function: contSpec
#
# Using a simplified L-curve method to compute the continuous relaxation 
# spectra H(s) given G*(w) from an input file
#
# Uses parameters from SetParameters function:
# These parameters are prefixed by par.XXXXX, where XXXXX is the parameter
# Definitions of these parameters are provided in the SetParameters function
#
# Output: Yields the fitted H(s), and optionally critical lambda "lamC"
# If verbose is on: rho-eta.dat (used to determine lamC)
#                   H.dat (the spectrum)
#                   Gfit.dat (G* implied by the extracted H(s))
#

contSpec <- function(par = NULL) {
  if (is.null(par)) {
    par <- SetParameters()  # Load in global settings
  }

  if (par$verbose) {
    cat("\n(*) Start\n(*) Loading Data File:", par$GstFile, "...")
  }

  exp_data <- GetExpData(par$GstFile)
  t <- exp_data$t
  Gexp <- exp_data$Gt

  if (par$verbose) {
    cat("done\n(*) Initial Set up...")
  }

  start_time <- Sys.time()

  n <- length(t)
  ns <- par$ns  # discretization of 'tau'

  tmin <- t[1]
  tmax <- max(t)

  smin <- 0
  smax <- 0

  switch(par.FreqEnd,
         `1` = {
           smin <- exp(-pi/2) * tmin
           smax <- exp(pi/2) * tmax
         },
         `2` = {
           smin <- tmin
           smax <- tmax
         },
         `3` = {
           smin <- exp(pi/2) * tmin
           smax <- exp(-pi/2) * tmax
         }
        ) 
  hs <- (smax / smin)^(1 / (ns - 1))
  s <- smin * hs^seq(0, ns - 1)

  Hgs <- InitializeH(Gexp, t, s)

  if (par$verbose) {
    cat("\n(*) Finding optimum lambda with 'lcurve'...")
  }

  if (par$lamC == 0) {
    lcurve_result <- lcurve(Gexp, Hgs, t, s, par$SmFacLam)
    lamC <- lcurve_result$lamC
    lam <- lcurve_result$lam
    rho <- lcurve_result$rho
    eta <- lcurve_result$eta
  } else {
    lamC <- par$lamC
  }

  if (par$verbose) {
    cat(sprintf("%e (%5.1f seconds)\n(*) Extracting the continuous spectrum, ...", lamC, as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
  }

  H <- LevenMarq(lamC, Gexp, Hgs, t, s)

  if (par$verbose) {
    cat(sprintf("done (%5.1f seconds)\n(*) Writing and Printing, ...", as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
  }

   
  if (par$verbose) {
    
    if (!dir.exists("output")) {
      if (!dir.create("output", showWarnings = FALSE)) {
        stop("Error: unable to create the 'output' directory.")
      }
    }
    
    if (par$lamC == 0) {
      write.table(data.frame(lam = lam, rho = rho, eta = eta), file = "output/rho-eta.dat", row.names = FALSE, col.names = FALSE, quote = FALSE)
    }

    write.table(data.frame(s = s, H = H), file = "output/H.dat", row.names = FALSE, col.names = FALSE, quote = FALSE)

    K <- kernel(H, t, s)
    Gfit <- data.frame(t = t, Gp_fit = K[1:n], Gpp_fit = K[(n + 1):(2 * n)])
    write.table(Gfit, file = "output/Gfit.dat", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }

  if (par$plotting) {
   
    plot(s, H, type = "o", log = "x", xlab = "s", ylab = "H(s)", main = "H")
    K <- kernel(H, t, s)
    plot(t, Gexp[1:n], type = "p", log = "xy", xlab = "t", ylab = "G*(exp), G*(fit)", main = "Gp", col = "red")
    points(t, K[1:n], type = "l", col = "black")
    points(t, Gexp[(n + 1):(2 * n)], type = "p", col = "blue")
    points(t, K[(n + 1):(2 * n)], type = "l", col = "black")
  }

  if (par$verbose) {
    cat("done\n(*) End\n")
  }

  return(list(H = H, lamC = lamC))
}
