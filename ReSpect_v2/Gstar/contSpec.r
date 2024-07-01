%# Function: contSpec
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
  w <- exp_data$w
  Gexp <- exp_data$Gexp

  if (par$verbose) {
    cat("done\n(*) Initial Set up...")
  }

  start_time <- Sys.time()

  n <- length(w)
  ns <- par$ns  # discretization of 'tau'

  wmin <- min(w)
  wmax <- max(w)

  smin <- 0
  smax <- 0

  switch(par$FreqEnd,
         `1` = {
           smin <- exp(-pi / 2) / wmax
           smax <- exp(pi / 2) / wmin
         },
         `2` = {
           smin <- 1 / wmax
           smax <- 1 / wmin
         },
         `3` = {
           smin <- exp(pi / 2) / wmax
           smax <- exp(-pi / 2) / wmin
         })

  hs <- (smax / smin)^(1 / (ns - 1))
  s <- smin * hs^seq(0, ns - 1)

  Hgs <- InitializeH(Gexp, w, s)

  if (par$verbose) {
    cat("\n(*) Finding optimum lambda with 'lcurve'...")
  }

  if (par$lamC == 0) {
    lcurve_result <- lcurve(Gexp, Hgs, w, s, par$SmFacLam)
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

  H <- LevenMarq(lamC, Gexp, Hgs, w, s)

  if (par$verbose) {
    cat(sprintf("done (%5.1f seconds)\n(*) Writing and Printing, ...", as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
  }

  if (par$verbose) {
    if (par$lamC == 0) {
      write.table(data.frame(lam = lam, rho = rho, eta = eta), file = "output/rho-eta.dat", row.names = FALSE, col.names = FALSE, quote = FALSE)
    }

    write.table(data.frame(s = s, H = H), file = "output/H.dat", row.names = FALSE, col.names = FALSE, quote = FALSE)

    K <- kernel(H, w, s)
    Gfit <- data.frame(w = w, Gp_fit = K[1:n], Gpp_fit = K[(n + 1):(2 * n)])
    write.table(Gfit, file = "output/Gfit.dat", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }

  if (par$plotting) {
    par(mfrow = c(2, 1))

    plot(s, H, type = "o", log = "x", xlab = "s", ylab = "H(s)", main = "H")
    K <- kernel(H, w, s)
    plot(w, Gexp[1:n], type = "p", log = "xy", xlab = "w", ylab = "G*(exp), G*(fit)", main = "Gp", col = "red")
    points(w, K[1:n], type = "l", col = "black")
    points(w, Gexp[(n + 1):(2 * n)], type = "p", col = "blue")
    points(w, K[(n + 1):(2 * n)], type = "l", col = "black")
  }

  if (par$verbose) {
    cat("done\n(*) End\n")
  }

  return(list(H = H, lamC = lamC))
}
