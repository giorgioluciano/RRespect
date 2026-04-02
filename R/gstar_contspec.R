# gstar_contspec.R — continuous relaxation spectrum, frequency domain (internal)
# Ported from Gstar/contSpec.R; no library() / source() calls.
# Uses .gstarGetExpData / .gstarGetKernMat / .gstarKernelPrestore from gstar_common.R.

# ── Jacobian helper ───────────────────────────────────────────────────────────
.gstarKernelD <- function(H, kernMat) {
  n  <- nrow(kernMat) / 2
  ns <- ncol(kernMat)
  Hsuper <- matrix(rep(exp(H), each = 2 * n), nrow = 2 * n, ncol = ns)
  kernMat * Hsuper
}

# ── LM residual (augmented with smoothness penalty) ───────────────────────────
.gstarResidualLM <- function(H, lam, Gexp, wexp, kernMat) {
  n  <- nrow(kernMat) / 2
  ns <- ncol(kernMat)
  nl <- ns - 2
  r  <- numeric(2 * n + nl)
  G0 <- NULL
  if (length(H) > ns) {
    G0         <- H[ns + 1]
    H          <- H[1:ns]
    r[1:(2*n)] <- wexp * (1 - .gstarKernelPrestore(H, kernMat, G0) / Gexp)
  } else {
    r[1:(2*n)] <- wexp * (1 - .gstarKernelPrestore(H, kernMat) / Gexp)
  }
  r[(2*n+1):(2*n+nl)] <- sqrt(lam) * diff(H, differences = 2)
  r
}

# ── LM Jacobian ───────────────────────────────────────────────────────────────
.gstarJacobianLM <- function(H, lam, Gexp, wexp, kernMat) {
  n  <- nrow(kernMat) / 2
  ns <- ncol(kernMat)
  nl <- ns - 2

  L <- -2 * diag(ns)
  L[row(L) == col(L) - 1] <- 1
  L[row(L) == col(L) + 1] <- 1
  L <- L[2:(nl + 1), ]

  Kmatrix <- (wexp / Gexp) %o% rep(1, ns)

  if (length(H) > ns) {
    G0 <- H[ns + 1]
    H  <- H[1:ns]
    Jr <- matrix(0, nrow = 2 * n + nl, ncol = ns + 1)
    Jr[1:(2*n),          1:ns]   <- -.gstarKernelD(H, kernMat) * Kmatrix
    Jr[1:n,              ns + 1] <- -wexp[1:n] / Gexp[1:n]
    Jr[(2*n+1):(2*n+nl), 1:ns]   <- sqrt(lam) * L
  } else {
    Jr <- matrix(0, nrow = 2 * n + nl, ncol = ns)
    Jr[1:(2*n),          1:ns] <- -.gstarKernelD(H, kernMat) * Kmatrix
    Jr[(2*n+1):(2*n+nl), 1:ns] <- sqrt(lam) * L
  }
  Jr
}

# ── getH: LM solve for H (±G0 plateau) ───────────────────────────────────────
.gstarGetH <- function(lam, Gexp, wexp, H, kernMat, G0 = NULL) {
  ctrl <- minpack.lm::nls.lm.control(maxiter = 500, ftol = 1e-8, ptol = 1e-8)
  if (!is.null(G0)) {
    Hplus <- c(H, G0)
    res   <- minpack.lm::nls.lm(
      par     = Hplus,
      fn      = .gstarResidualLM,
      jac     = .gstarJacobianLM,
      lam     = lam, Gexp = Gexp, wexp = wexp, kernMat = kernMat,
      control = ctrl
    )
    return(list(H = res$par[seq_len(length(H))],
                G0 = res$par[length(H) + 1L]))
  }
  res <- minpack.lm::nls.lm(
    par     = H,
    fn      = .gstarResidualLM,
    jac     = .gstarJacobianLM,
    lam     = lam, Gexp = Gexp, wexp = wexp, kernMat = kernMat,
    control = ctrl
  )
  res$par
}

# ── Regularisation matrices ───────────────────────────────────────────────────
.gstarGetAmatrix <- function(ns) {
  nl  <- ns - 2
  idx <- seq_len(nl)
  L   <- matrix(0, nl, ns)
  L[cbind(idx, idx    )] <-  1
  L[cbind(idx, idx + 1)] <- -2
  L[cbind(idx, idx + 2)] <-  1
  t(L) %*% L
}

.gstarGetBmatrix <- function(H, kernMat, Gexp, wexp, G0 = NULL) {
  ns      <- length(H)
  Kmatrix <- (wexp / Gexp) %o% rep(1, ns)
  Jr      <- -.gstarKernelD(H, kernMat) * Kmatrix
  r       <- if (!is.null(G0)) {
    wexp * (1 - .gstarKernelPrestore(H, kernMat, G0) / Gexp)
  } else {
    wexp * (1 - .gstarKernelPrestore(H, kernMat) / Gexp)
  }
  t(Jr) %*% Jr + diag(as.vector(t(r) %*% Jr))
}

# ── InitializeH ───────────────────────────────────────────────────────────────
.gstarInitializeH <- function(Gexp, wexp, s, kernMat, G0 = NULL) {
  H   <- -5.0 * rep(1, length(s)) + sin(pi * s)
  lam <- 1e0
  if (!is.null(G0)) {
    return(.gstarGetH(lam, Gexp, wexp, H, kernMat, G0))
  }
  .gstarGetH(lam, Gexp, wexp, H, kernMat)
}

# ── L-curve (Bayesian, frequency domain) ─────────────────────────────────────
.gstarLcurve <- function(Gexp, wexp, Hgs, kernMat, par, G0 = NULL) {
  npoints <- as.integer(par$lamDensity * (log10(par$lam_max) - log10(par$lam_min)))
  hlam    <- (par$lam_max / par$lam_min)^(1 / (npoints - 1))
  lam     <- par$lam_min * hlam^(0:(npoints - 1))

  ns      <- length(Hgs)
  eta     <- numeric(npoints)
  rho     <- numeric(npoints)
  logP    <- numeric(npoints)
  Hlambda <- matrix(0, nrow = ns, ncol = npoints)
  Amat    <- .gstarGetAmatrix(ns)
  LogDetN <- as.numeric(determinant(Amat, logarithm = TRUE)$modulus)

  H       <- Hgs
  logPmax <- -Inf
  i_start <- 1

  for (i in rev(seq_len(npoints))) {
    lamb <- lam[i]
    if (isTRUE(par$plateau)) {
      res    <- .gstarGetH(lamb, Gexp, wexp, H, kernMat, G0)
      H <- res$H;  G0 <- res$G0
      rho[i] <- sqrt(sum((wexp * (1 - .gstarKernelPrestore(H, kernMat, G0) / Gexp))^2))
      Bmat   <- .gstarGetBmatrix(H, kernMat, Gexp, wexp, G0)
    } else {
      H      <- .gstarGetH(lamb, Gexp, wexp, H, kernMat)
      rho[i] <- sqrt(sum((wexp * (1 - .gstarKernelPrestore(H, kernMat) / Gexp))^2))
      Bmat   <- .gstarGetBmatrix(H, kernMat, Gexp, wexp)
    }
    eta[i]       <- sqrt(sum(diff(H, differences = 2)^2))
    Hlambda[, i] <- H

    LogDetC <- as.numeric(determinant(lamb * Amat + Bmat, logarithm = TRUE)$modulus)
    V       <- rho[i]^2 + lamb * eta[i]^2
    logP[i] <- -V + 0.5 * (LogDetN + ns * log(lamb) - LogDetC) - lamb

    if (logP[i] > logPmax) {
      logPmax <- logP[i]
    } else if (logP[i] < logPmax - 18) {
      i_start <- i
      break
    }
  }

  idx     <- i_start:npoints
  lam     <- lam[idx];  logP <- logP[idx]
  eta     <- eta[idx];  rho  <- rho[idx]
  Hlambda <- Hlambda[, idx, drop = FALSE]
  logP    <- logP - max(logP)

  plam <- exp(logP);  plam <- plam / sum(plam)
  lamM <- exp(sum(plam * log(lam)))

  if (par$SmFacLam > 0) {
    lamM <- exp(log(lamM) + par$SmFacLam * (max(log(lam)) - log(lamM)))
  } else if (par$SmFacLam < 0) {
    lamM <- exp(log(lamM) + par$SmFacLam * (log(lamM) - min(log(lam))))
  }

  if (par$verbose) cat(sprintf("L-curve: lamM = %.4e\n", lamM))

  list(lamM = lamM, lam = lam, rho = rho, eta = eta,
       logP = logP, Hlambda = Hlambda, G0 = G0)
}

# ── contSpec (main entry point) ───────────────────────────────────────────────
.gstarContSpec <- function(par, writeOutput = FALSE, outputDir = "output") {
  if (par$verbose)
    cat(sprintf("\n(*) Start\n(*) Loading Data File: %s...\n", par$GstFile))

  inp  <- .gstarGetExpData(par$GstFile)
  w    <- inp$w;  Gexp <- inp$Gst;  wexp <- inp$wexp
  n    <- length(w)
  ns   <- par$ns

  wmin <- w[1];  wmax <- w[n]

  if (par$FreqEnd == 1) {
    smin <- exp(-pi / 2) / wmax;  smax <- exp(pi / 2)  / wmin
  } else if (par$FreqEnd == 2) {
    smin <- 1 / wmax;              smax <- 1 / wmin
  } else {
    smin <- exp(pi / 2) / wmax;   smax <- exp(-pi / 2) / wmin
  }

  hs      <- (smax / smin)^(1 / (ns - 1))
  s       <- smin * hs^(0:(ns - 1))
  kernMat <- .gstarGetKernMat(s, w)

  if (par$verbose) cat("(*) Initial Set up...\n")

  if (isTRUE(par$plateau)) {
    res <- .gstarInitializeH(Gexp, wexp, s, kernMat, min(Gexp))
    Hgs <- res$H;  G0 <- res$G0
  } else {
    Hgs <- .gstarInitializeH(Gexp, wexp, s, kernMat)
    G0  <- NULL
  }

  if (par$verbose) cat("(*) Building the L-curve ...\n")

  if (par$lamC == 0) {
    lc   <- .gstarLcurve(Gexp, wexp, Hgs, kernMat, par, G0)
    lamC <- lc$lamM
    lam  <- lc$lam;  rho <- lc$rho;  eta <- lc$eta
    logP <- lc$logP;  Hlam <- lc$Hlambda
    if (isTRUE(par$plateau)) G0 <- lc$G0
  } else {
    lamC <- par$lamC
  }

  if (par$verbose)
    cat(sprintf("(*) Extracting CRS ... lamM = %.3e\n", lamC))

  if (isTRUE(par$plateau)) {
    res <- .gstarGetH(lamC, Gexp, wexp, Hgs, kernMat, G0)
    H   <- res$H;  G0 <- res$G0
    cat(sprintf("G0 = %.3e\n", G0))
  } else {
    H  <- .gstarGetH(lamC, Gexp, wexp, Hgs, kernMat)
    G0 <- 0
  }

  if (isTRUE(writeOutput)) {
    if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)
    if (isTRUE(par$plateau)) {
      K <- .gstarKernelPrestore(H, kernMat, G0)
      write(sprintf("# G0 = %e", G0), file = file.path(outputDir, "H.dat"))
      utils::write.table(cbind(s, H, rep(0, ns)), file.path(outputDir, "H.dat"),
                         append = TRUE, row.names = FALSE, col.names = FALSE)
    } else {
      K <- .gstarKernelPrestore(H, kernMat)
      utils::write.table(cbind(s, H, rep(0, ns)), file.path(outputDir, "H.dat"),
                         row.names = FALSE, col.names = FALSE)
    }
    utils::write.table(cbind(w, K[1:n], K[(n+1):(2*n)]),
                       file.path(outputDir, "Gfit.dat"),
                       row.names = FALSE, col.names = FALSE)

    if (par$lamC == 0) {
      fhlam <- file(file.path(outputDir, "Hlam.dat"), "wt")
      for (i in seq_along(lam))
        writeLines(paste(Hlam[, i], collapse = " "), fhlam)
      close(fhlam)
      utils::write.table(cbind(lam, logP), file.path(outputDir, "logPlam.dat"),
                         row.names = FALSE, col.names = FALSE)
      utils::write.table(cbind(lam, rho, eta), file.path(outputDir, "rho-eta.dat"),
                         row.names = FALSE, col.names = FALSE)
    }
  }

  if (par$verbose) cat("(*) Done\n")

  list(s = s, H = H, G0 = G0, lamC = lamC,
       w = w, Gexp = Gexp, wexp = wexp, kernMat = kernMat)
}
