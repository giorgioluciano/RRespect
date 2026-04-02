# gt_contspec.R — continuous relaxation spectrum, time domain (internal)
# Ported from Gt/contSpec.R; no library() / source() calls.

.gtInitializeH <- function(Gexp, wexp, s, kernMat) {
  H   <- -5 * rep(1, length(s)) + sin(pi * s)
  lam <- 1e0
  .gtLevenMarq(lam, Gexp, wexp, H, kernMat)
}

.gtContSpec <- function(par, writeOutput = FALSE, outputDir = "output") {
  if (par$verbose) cat("\n(*) Start\n(*) Loading Data File:", par$GtFile, "\n")

  dat  <- .gtGetExpData(par$GtFile)
  t    <- dat$t;  Gexp <- dat$Gt;  wexp <- dat$wexp
  n    <- length(t)
  ns   <- par$ns

  tmin <- t[1];  tmax <- t[n]

  if (par$FreqEnd == 1) {
    smin <- exp(-pi / 2) * tmin;  smax <- exp(pi / 2) * tmax
  } else if (par$FreqEnd == 2) {
    smin <- tmin;  smax <- tmax
  } else {
    smin <- exp(pi / 2) * tmin;  smax <- exp(-pi / 2) * tmax
  }

  hs <- (smax / smin)^(1 / (ns - 1))
  s  <- smin * hs^(0:(ns - 1))

  kernMat <- .gtGetKernMat(s, t)

  if (par$verbose) cat("(*) Initializing H...\n")
  Hgs <- .gtInitializeH(Gexp, wexp, s, kernMat)

  if (par$lamC == 0) {
    if (par$verbose) cat("(*) Building L-curve (Bayesian)...\n")
    res  <- .gtLcurve(Gexp, wexp, Hgs, kernMat, par)
    lamC <- res$lamC
    lam  <- res$lam;  rho <- res$rho;  eta <- res$eta
    logP <- res$logP
  } else {
    lamC <- par$lamC
  }

  if (par$verbose) cat(sprintf("(*) lamM = %.3e\n", lamC))

  H    <- .gtLevenMarq(lamC, Gexp, wexp, Hgs, kernMat)
  Kfit <- .gtKernelPrestore(H, kernMat)

  if (isTRUE(writeOutput)) {
    if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)
    utils::write.table(cbind(s, H),     file.path(outputDir, "H.dat"),
                       row.names = FALSE, col.names = FALSE)
    utils::write.table(cbind(t, Kfit),  file.path(outputDir, "Gfit.dat"),
                       row.names = FALSE, col.names = FALSE)
    if (par$lamC == 0) {
      utils::write.table(cbind(lam, rho, eta), file.path(outputDir, "rho-eta.dat"),
                         row.names = FALSE, col.names = FALSE)
    }
  }

  if (par$verbose) cat(sprintf("(*) CRS done. lamC = %.3e\n", lamC))

  if (par$plotting && isTRUE(writeOutput)) {
    grDevices::pdf(file.path(outputDir, "H_time.pdf"))
    graphics::plot(s, H, type = "l", log = "x",
                   xlab = expression(tau), ylab = "H(tau)",
                   main = "Continuous Relaxation Spectrum (time)")
    grDevices::dev.off()

    grDevices::pdf(file.path(outputDir, "Gfit_time.pdf"))
    graphics::plot(t, Gexp, log = "xy", col = "blue", pch = 19, cex = 0.5,
                   xlab = "t", ylab = "G(t)")
    graphics::lines(t, Kfit, col = "black", lwd = 2)
    grDevices::dev.off()
  }

  list(s = s, H = H, lamC = lamC,
       kernMat = kernMat, t = t, Gexp = Gexp, wexp = wexp)
}
