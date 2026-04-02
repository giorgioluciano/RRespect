# gt_lcurve.R — Bayesian L-curve regularisation, time domain (internal)
# Ported from Gt/LCurve.R; no library() / source() calls.

.gtGetAmatrix <- function(ns) {
  nl  <- ns - 2
  idx <- seq_len(nl)
  L   <- matrix(0, nl, ns)
  L[cbind(idx, idx    )] <-  1
  L[cbind(idx, idx + 1)] <- -2
  L[cbind(idx, idx + 2)] <-  1
  t(L) %*% L
}

.gtGetBmatrix <- function(H, kernMat, Gexp, wexp) {
  n  <- nrow(kernMat);  ns <- ncol(kernMat)
  Kmatrix <- matrix(wexp / Gexp, nrow = n, ncol = ns)
  Jr <- -.gtKernelD(H, kernMat) * Kmatrix
  r  <- wexp * (1 - .gtKernelPrestore(H, kernMat) / Gexp)
  t(Jr) %*% Jr + diag(as.vector(t(r) %*% Jr))
}

.gtLcurve <- function(Gexp, wexp, Hgs, kernMat, par) {
  ns      <- ncol(kernMat)
  npoints <- as.integer(par$lamDensity * (log10(par$lam_max) - log10(par$lam_min)))
  hlam    <- (par$lam_max / par$lam_min)^(1 / (npoints - 1))
  lam     <- par$lam_min * hlam^(0:(npoints - 1))

  eta  <- numeric(npoints);  rho  <- numeric(npoints)
  logP <- numeric(npoints);  Hlambda <- matrix(0, ns, npoints)

  Amat    <- .gtGetAmatrix(ns)
  LogDetN <- as.numeric(determinant(Amat, logarithm = TRUE)$modulus)

  H       <- Hgs
  logPmax <- -Inf
  i_break <- 1

  for (i in rev(seq_len(npoints))) {
    lamb <- lam[i]
    H    <- .gtLevenMarq(lamb, Gexp, wexp, H, kernMat)

    rho[i]  <- sqrt(sum((wexp * (1 - .gtKernelPrestore(H, kernMat) / Gexp))^2))
    eta[i]  <- sqrt(sum(diff(H, differences = 2)^2))
    Hlambda[, i] <- H

    Bmat    <- .gtGetBmatrix(H, kernMat, Gexp, wexp)
    LogDetC <- as.numeric(determinant(lamb * Amat + Bmat, logarithm = TRUE)$modulus)
    V       <- rho[i]^2 + lamb * eta[i]^2
    logP[i] <- -V + 0.5 * (LogDetN + ns * log(lamb) - LogDetC) - lamb

    if (logP[i] > logPmax) {
      logPmax <- logP[i]
    } else if (logP[i] < logPmax - 18) {
      i_break <- i
      break
    }
  }

  idx <- i_break:npoints
  lam     <- lam[idx];  logP <- logP[idx]
  rho     <- rho[idx];  eta  <- eta[idx]
  Hlambda <- Hlambda[, idx, drop = FALSE]
  logP    <- logP - max(logP)

  plam <- exp(logP);  plam <- plam / sum(plam)
  lamM <- exp(sum(plam * log(lam)))

  SmFac <- par$SmFacLam
  if (SmFac > 0)      lamM <- exp(log(lamM) + SmFac * (max(log(lam)) - log(lamM)))
  else if (SmFac < 0) lamM <- exp(log(lamM) + SmFac * (log(lamM) - min(log(lam))))

  list(lamC = lamM, lam = lam, rho = rho, eta = eta,
       logP = logP, Hlambda = Hlambda)
}
