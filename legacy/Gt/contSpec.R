
# contSpec_time.R
# Continuous Relaxation Spectrum — time domain G(t)
# Faithfully translated from contSpec.py (pyReSpect-time)

source("LCurve.R")

InitializeH <- function(Gexp, wexp, s, kernMat) {
  H   <- -5 * rep(1, length(s)) + sin(pi * s)
  lam <- 1e0
  H   <- LevenMarq(lam, Gexp, wexp, H, kernMat)
  return(H)
}

contSpec <- function(par) {
  if (missing(par)) par <- SetParameters()

  if (par$verbose) cat("\n(*) Start\n(*) Loading Data File:", par$GtFile, "\n")

  dat  <- GetExpData(par$GtFile)
  t    <- dat$t;  Gexp <- dat$Gt;  wexp <- dat$wexp
  n    <- length(t)
  ns   <- par$ns

  tmin <- t[1];  tmax <- t[n]

  if (par$FreqEnd == 1) {
    smin <- exp(-pi/2) * tmin;  smax <- exp(pi/2) * tmax
  } else if (par$FreqEnd == 2) {
    smin <- tmin;  smax <- tmax
  } else {
    smin <- exp(pi/2) * tmin;  smax <- exp(-pi/2) * tmax
  }

  hs <- (smax / smin)^(1 / (ns - 1))
  s  <- smin * hs^(0:(ns-1))

  kernMat <- getKernMat(s, t)

  if (par$verbose) cat("(*) Initializing H...\n")
  Hgs <- InitializeH(Gexp, wexp, s, kernMat)

  if (par$lamC == 0) {
    if (par$verbose) cat("(*) Building L-curve (Bayesian)...\n")
    res  <- lcurve(Gexp, wexp, Hgs, kernMat, par)
    lamC <- res$lamC
    lam  <- res$lam;  rho <- res$rho;  eta <- res$eta
    logP <- res$logP
  } else {
    lamC <- par$lamC
  }

  if (par$verbose) cat(sprintf("(*) lamM = %.3e\n", lamC))

  H <- LevenMarq(lamC, Gexp, wexp, Hgs, kernMat)

  if (!dir.exists("output")) dir.create("output")
  Kfit <- kernel_prestore(H, kernMat)

  write.table(cbind(s, H), "output/H.dat",
              row.names = FALSE, col.names = FALSE)
  write.table(cbind(t, Kfit), "output/Gfit.dat",
              row.names = FALSE, col.names = FALSE)

  if (par$lamC == 0) {
    write.table(cbind(lam, rho, eta), "output/rho-eta.dat",
                row.names = FALSE, col.names = FALSE)
  }

  if (par$verbose) cat(sprintf("(*) CRS done. lamC = %.3e\n", lamC))

  if (par$plotting) {
    pdf("output/H_time.pdf")
    plot(s, H, type = "l", log = "x",
         xlab = expression(tau), ylab = "H(tau)",
         main = "Continuous Relaxation Spectrum (time)")
    dev.off()

    pdf("output/Gfit_time.pdf")
    plot(t, Gexp, log = "xy", col = "blue", pch = 19, cex = 0.5,
         xlab = "t", ylab = "G(t)")
    lines(t, Kfit, col = "black", lwd = 2)
    dev.off()
  }

  return(list(s = s, H = H, lamC = lamC,
              kernMat = kernMat, t = t, Gexp = Gexp, wexp = wexp))
}
