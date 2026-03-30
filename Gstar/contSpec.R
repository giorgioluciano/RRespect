# contSpec.R
# Traduzione fedele di contSpec.py (pyReSpect-freq)
# Firma: contSpec(par) -- legge il .dat da par$GstFile internamente

# ── Helper: InitializeH ───────────────────────────────────────────────────────
InitializeH <- function(Gexp, wexp, s, kernMat, G0 = NULL) {
  H   <- -5.0 * rep(1, length(s)) + sin(pi * s)
  lam <- 1e0
  if (!is.null(G0)) {
    res <- getH(lam, Gexp, wexp, H, kernMat, G0)
    return(list(H = res$H, G0 = res$G0))
  } else {
    H <- getH(lam, Gexp, wexp, H, kernMat)
    return(H)
  }
}

# ── Helper: getAmatrix ────────────────────────────────────────────────────────
getAmatrix <- function(ns) {
  nl <- ns - 2
  L  <- -2 * diag(ns)
  L[row(L) == col(L) - 1] <- 1
  L[row(L) == col(L) + 1] <- 1
  L  <- L[2:(nl+1), ]
  return(t(L) %*% L)
}

# ── Helper: getBmatrix ────────────────────────────────────────────────────────
getBmatrix <- function(H, kernMat, Gexp, wexp, G0 = NULL) {
  n  <- nrow(kernMat) / 2
  ns <- length(H)
  Kmatrix <- (wexp / Gexp) %o% rep(1, ns)
  Jr <- -kernelD(H, kernMat) * Kmatrix
  if (!is.null(G0)) {
    r <- wexp * (1 - kernel_prestore(H, kernMat, G0) / Gexp)
  } else {
    r <- wexp * (1 - kernel_prestore(H, kernMat) / Gexp)
  }
  B <- t(Jr) %*% Jr + diag(as.vector(t(r) %*% Jr))
  return(B)
}

# ── Helper: kernelD ───────────────────────────────────────────────────────────
kernelD <- function(H, kernMat) {
  n      <- nrow(kernMat) / 2
  ns     <- ncol(kernMat)
  Hsuper <- matrix(rep(exp(H), each = 2*n), nrow = 2*n, ncol = ns)
  return(kernMat * Hsuper)
}

# ── Helper: residualLM ────────────────────────────────────────────────────────
residualLM <- function(H, lam, Gexp, wexp, kernMat) {
  n  <- nrow(kernMat) / 2
  ns <- ncol(kernMat)
  nl <- ns - 2
  r  <- numeric(2*n + nl)
  G0 <- NULL
  if (length(H) > ns) {
    G0   <- H[ns + 1]
    H    <- H[1:ns]
    r[1:(2*n)] <- wexp * (1 - kernel_prestore(H, kernMat, G0) / Gexp)
  } else {
    r[1:(2*n)] <- wexp * (1 - kernel_prestore(H, kernMat) / Gexp)
  }
  r[(2*n+1):(2*n+nl)] <- sqrt(lam) * diff(H, differences = 2)
  return(r)
}

# ── Helper: jacobianLM ────────────────────────────────────────────────────────
jacobianLM <- function(H, lam, Gexp, wexp, kernMat) {
  n  <- nrow(kernMat) / 2
  ns <- ncol(kernMat)
  nl <- ns - 2
  L  <- -2 * diag(ns)
  L[row(L) == col(L) - 1] <- 1
  L[row(L) == col(L) + 1] <- 1
  L  <- L[2:(nl+1), ]
  Kmatrix <- (wexp / Gexp) %o% rep(1, ns)
  if (length(H) > ns) {
    G0 <- H[ns + 1]
    H  <- H[1:ns]
    Jr <- matrix(0, nrow = 2*n + nl, ncol = ns + 1)
    Jr[1:(2*n),  1:ns]   <- -kernelD(H, kernMat) * Kmatrix
    Jr[1:n,      ns+1]   <- -wexp[1:n] / Gexp[1:n]
    Jr[(2*n+1):(2*n+nl), 1:ns] <- sqrt(lam) * L
  } else {
    Jr <- matrix(0, nrow = 2*n + nl, ncol = ns)
    Jr[1:(2*n),  1:ns]   <- -kernelD(H, kernMat) * Kmatrix
    Jr[(2*n+1):(2*n+nl), 1:ns] <- sqrt(lam) * L
  }
  return(Jr)
}

# ── Helper: getH ─────────────────────────────────────────────────────────────
getH <- function(lam, Gexp, wexp, H, kernMat, G0 = NULL) {
  if (!is.null(G0)) {
    Hplus <- c(H, G0)
    res   <- nls.lm(
      par     = Hplus,
      fn      = residualLM,
      jac     = jacobianLM,
      lam     = lam, Gexp = Gexp, wexp = wexp, kernMat = kernMat,
      control = nls.lm.control(maxiter = 500, ftol = 1e-8, ptol = 1e-8)
    )
    return(list(H = res$par[1:length(H)], G0 = res$par[length(H)+1]))
  } else {
    res <- nls.lm(
      par     = H,
      fn      = residualLM,
      jac     = jacobianLM,
      lam     = lam, Gexp = Gexp, wexp = wexp, kernMat = kernMat,
      control = nls.lm.control(maxiter = 500, ftol = 1e-8, ptol = 1e-8)
    )
    return(res$par)
  }
}

# ── lcurve ────────────────────────────────────────────────────────────────────
lcurve <- function(Gexp, wexp, Hgs, kernMat, par, G0 = NULL) {
  library(Matrix)
  npoints <- as.integer(par$lamDensity * (log10(par$lam_max) - log10(par$lam_min)))
  hlam    <- (par$lam_max / par$lam_min)^(1 / (npoints - 1))
  lam     <- par$lam_min * hlam^(0:(npoints-1))
  eta     <- numeric(npoints)
  rho     <- numeric(npoints)
  logP    <- numeric(npoints)
  H       <- Hgs
  ns      <- length(H)
  Hlambda <- matrix(0, nrow = ns, ncol = npoints)
  Amat    <- getAmatrix(ns)
  LogDetN <- as.numeric(determinant(Amat, logarithm = TRUE)$modulus)
  logPmax <- -Inf
  i_start <- 1

  for (i in rev(seq_along(lam))) {
    lamb <- lam[i]
    if (isTRUE(par$plateau)) {
      res    <- getH(lamb, Gexp, wexp, H, kernMat, G0)
      H <- res$H; G0 <- res$G0
      rho[i] <- sqrt(sum((wexp * (1 - kernel_prestore(H, kernMat, G0) / Gexp))^2))
      Bmat   <- getBmatrix(H, kernMat, Gexp, wexp, G0)
    } else {
      H      <- getH(lamb, Gexp, wexp, H, kernMat)
      rho[i] <- sqrt(sum((wexp * (1 - kernel_prestore(H, kernMat) / Gexp))^2))
      Bmat   <- getBmatrix(H, kernMat, Gexp, wexp)
    }
    eta[i]        <- sqrt(sum(diff(H, differences = 2)^2))
    Hlambda[, i]  <- H
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

  lam     <- lam[i_start:npoints]
  logP    <- logP[i_start:npoints]
  eta     <- eta[i_start:npoints]
  rho     <- rho[i_start:npoints]
  Hlambda <- Hlambda[, i_start:npoints, drop = FALSE]
  logP    <- logP - max(logP)

  plam <- exp(logP); plam <- plam / sum(plam)
  lamM <- exp(sum(plam * log(lam)))

  if (par$SmFacLam > 0) {
    lamM <- exp(log(lamM) + par$SmFacLam * (max(log(lam)) - log(lamM)))
  } else if (par$SmFacLam < 0) {
    lamM <- exp(log(lamM) + par$SmFacLam * (log(lamM) - min(log(lam))))
  }

  if (par$verbose)
    cat(sprintf("L-curve: lamM = %.4e\n", lamM))

  return(list(lamM = lamM, lam = lam, rho = rho, eta = eta,
              logP = logP, Hlambda = Hlambda, G0 = G0))
}

# ── contSpec (traduzione di getContSpec) ──────────────────────────────────────
contSpec <- function(par) {
  library(minpack.lm)

  if (par$verbose)
    cat(sprintf("\n(*) Start\n(*) Loading Data File: %s...\n", par$GstFile))

  inp     <- GetExpData(par$GstFile)
  w       <- inp$w
  Gexp    <- inp$Gst
  wexp    <- inp$wexp
  n       <- length(w)
  ns      <- par$ns

  wmin <- w[1]; wmax <- w[n]

  if (par$FreqEnd == 1) {
    smin <- exp(-pi/2) / wmax; smax <- exp(pi/2) / wmin
  } else if (par$FreqEnd == 2) {
    smin <- 1 / wmax;          smax <- 1 / wmin
  } else if (par$FreqEnd == 3) {
    smin <- exp(pi/2) / wmax;  smax <- exp(-pi/2) / wmin
  }

  hs      <- (smax / smin)^(1 / (ns - 1))
  s       <- smin * hs^(0:(ns-1))
  kernMat <- getKernMat(s, w)

  if (par$verbose) cat("(*) Initial Set up...\n")

  if (isTRUE(par$plateau)) {
    res  <- InitializeH(Gexp, wexp, s, kernMat, min(Gexp))
    Hgs  <- res$H; G0 <- res$G0
  } else {
    Hgs  <- InitializeH(Gexp, wexp, s, kernMat)
    G0   <- NULL
  }

  if (par$verbose) cat("(*) Building the L-curve ...\n")

  if (par$lamC == 0) {
    lc   <- lcurve(Gexp, wexp, Hgs, kernMat, par, G0)
    lamC <- lc$lamM
    lam  <- lc$lam; rho <- lc$rho; eta <- lc$eta
    logP <- lc$logP; Hlam <- lc$Hlambda
    if (isTRUE(par$plateau)) G0 <- lc$G0
  } else {
    lamC <- par$lamC
  }

  if (par$verbose)
    cat(sprintf("(*) Extracting CRS ... lamM = %.3e\n", lamC))

  if (isTRUE(par$plateau)) {
    res <- getH(lamC, Gexp, wexp, Hgs, kernMat, G0)
    H   <- res$H; G0 <- res$G0
    cat(sprintf("G0 = %.3e\n", G0))
  } else {
    H <- getH(lamC, Gexp, wexp, Hgs, kernMat)
    G0 <- 0
  }

  # Output files
  if (par$verbose) {
    if (isTRUE(par$plateau)) {
      K <- kernel_prestore(H, kernMat, G0)
      write(sprintf("# G0 = %e", G0), file = "output/H.dat")
      write.table(cbind(s, H, rep(0, ns)), "output/H.dat",
                  append = TRUE, row.names = FALSE, col.names = FALSE)
    } else {
      K <- kernel_prestore(H, kernMat)
      write.table(cbind(s, H, rep(0, ns)), "output/H.dat",
                  row.names = FALSE, col.names = FALSE)
    }
    write.table(cbind(w, K[1:n], K[(n+1):(2*n)]), "output/Gfit.dat",
                row.names = FALSE, col.names = FALSE)

    if (par$lamC == 0) {
      fhlam <- file("output/Hlam.dat", "wt")
      for (i in seq_along(lam))
        writeLines(paste(Hlam[, i], collapse = " "), fhlam)
      close(fhlam)
      write.table(cbind(lam, logP), "output/logPlam.dat",
                  row.names = FALSE, col.names = FALSE)
      write.table(cbind(lam, rho, eta), "output/rho-eta.dat",
                  row.names = FALSE, col.names = FALSE)
    }
    cat("(*) Done\n")
  }

  return(list(s = s, H = H, G0 = G0, lamC = lamC))
}
