
# LCurve_time.R  — Bayesian Hansen method, time-domain version
# Difference vs freq: kernMat is n x ns → getBmatrix uses n-dim residual

source("LevenMarq.R")

getAmatrix <- function(ns) {
  nl <- ns - 2
  L  <- diag(-2, ns)
  L[row(L) == col(L) - 1] <- 1
  L[row(L) == col(L) + 1] <- 1
  L_sub <- L[2:(nl+1), ]
  return(t(L_sub) %*% L_sub)
}

getBmatrix_time <- function(H, kernMat, Gexp, wexp) {
  n  <- nrow(kernMat);  ns <- ncol(kernMat)
  Kmatrix <- matrix(wexp / Gexp, nrow = n, ncol = ns)
  Jr <- -kernelD_time(H, kernMat) * Kmatrix
  r  <- wexp * (1 - kernel_prestore(H, kernMat) / Gexp)
  return(t(Jr) %*% Jr + diag(as.vector(t(r) %*% Jr)))
}

lcurve <- function(Gexp, wexp, Hgs, kernMat, par) {
  ns      <- ncol(kernMat)
  npoints <- as.integer(par$lamDensity * (log10(par$lam_max) - log10(par$lam_min)))
  hlam    <- (par$lam_max / par$lam_min)^(1 / (npoints - 1))
  lam     <- par$lam_min * hlam^(0:(npoints-1))

  eta  <- numeric(npoints);  rho  <- numeric(npoints)
  logP <- numeric(npoints);  Hlambda <- matrix(0, ns, npoints)

  Amat    <- getAmatrix(ns)
  LogDetN <- as.numeric(determinant(Amat, logarithm = TRUE)$modulus)

  H       <- Hgs
  logPmax <- -Inf
  
  i_break <- 1 
  
  for (i in rev(seq_len(npoints))) {
    lamb <- lam[i]
    H    <- LevenMarq(lamb, Gexp, wexp, H, kernMat)

    rho[i]  <- sqrt(sum((wexp * (1 - kernel_prestore(H, kernMat) / Gexp))^2))
    eta[i]  <- sqrt(sum(diff(H, differences = 2)^2))
    Hlambda[, i] <- H

    Bmat    <- getBmatrix_time(H, kernMat, Gexp, wexp)
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
  lam  <- lam[idx];  logP <- logP[idx]
  rho  <- rho[idx];  eta  <- eta[idx]
  Hlambda <- Hlambda[, idx, drop = FALSE]
  logP <- logP - max(logP)

  plam <- exp(logP);  plam <- plam / sum(plam)
  lamM <- exp(sum(plam * log(lam)))

  SmFac <- par$SmFacLam
  if (SmFac > 0)      lamM <- exp(log(lamM) + SmFac * (max(log(lam)) - log(lamM)))
  else if (SmFac < 0) lamM <- exp(log(lamM) + SmFac * (log(lamM) - min(log(lam))))

  return(list(lamC = lamM, lam = lam, rho = rho, eta = eta,
              logP = logP, Hlambda = Hlambda))
}
