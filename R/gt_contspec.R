# gt_contspec.R — continuous relaxation spectrum, time domain (internal)
# Ported from Gt/contSpec.R; no library() / source() calls.

.gtInitializeH <- function(
  Gexp, wexp, s, kernMat, G0 = NULL, solver = "legacy"
) {
  out <- .gtInitializeFit(Gexp, wexp, s, kernMat, G0, solver)
  if (is.null(G0)) out$H else list(H = out$H, G0 = out$G0)
}

.gtContSpec <- function(par, writeOutput = FALSE, outputDir = "output") {
  if (par$verbose) cat("\n(*) Start\n(*) Loading Data File:", par$GtFile, "\n")

  dat  <- .gtGetExpData(par$GtFile)
  t    <- dat$t;  Gexp <- dat$Gt;  wexp <- dat$wexp
  n    <- length(t)
  ns   <- par$ns
  solver <- if (is.null(par$solver)) "legacy" else par$solver
  scan_diagnostics <- NULL
  scan_result <- NULL

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
  initial <- .gtInitializeFit(
    Gexp, wexp, s, kernMat,
    G0 = if (isTRUE(par$plateau)) min(Gexp) else NULL,
    solver = solver
  )
  Hgs <- initial$H
  G0_initial <- initial$G0

  if (par$lamC == 0) {
    if (par$verbose) cat("(*) Building L-curve (Bayesian)...\n")
    res <- .gtLcurve(
      Gexp, wexp, Hgs, kernMat, par, G0_initial
    )
    lamC <- res$lamC
    lam  <- res$lam;  rho <- res$rho;  eta <- res$eta
    logP <- res$logP
    scan_diagnostics <- res$diagnostics
    scan_result <- res
  } else {
    lamC <- par$lamC
  }

  if (par$verbose) cat(sprintf("(*) lamM = %.3e\n", lamC))

  final <- .gtFitH(
    lamC, Gexp, wexp, Hgs, kernMat,
    G0 = G0_initial,
    solver = solver
  )
  H <- final$H
  G0 <- final$G0

  Kfit <- .gtKernelPrestore(H, kernMat, G0)

  if (isTRUE(writeOutput)) {
    if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)
	    h_file <- file.path(outputDir, "H.dat")

    if (isTRUE(par$plateau)) {
      write(sprintf("# G0 = %.17e", G0), h_file)
      utils::write.table(
        cbind(s, H), h_file,
        append = TRUE,
        row.names = FALSE,
        col.names = FALSE
      )
    } else {
      utils::write.table(
        cbind(s, H), h_file,
        row.names = FALSE,
        col.names = FALSE
      )
    }
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

    list(
    s = s,
    H = H,
    G0 = G0,
    lamC = lamC,
    kernMat = kernMat,
    t = t,
    Gexp = Gexp,
    wexp = wexp,
    solver = final$solver,
    scan_result = scan_result,
    diagnostics = list(
      initial = initial$diagnostics,
      scan = scan_diagnostics,
      final = final$diagnostics
    )
  )

}
