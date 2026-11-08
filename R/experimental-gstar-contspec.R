
# Experimental frequency-domain continuous workflow.
# Reference: shane5ul/pyReSpect-freq, Python 2.1.0,
# commit 908559720322941622ab39ced82e376776a39183.

.gstarGetHExperimental <- function(
  lam, Gexp, wexp, H, kernMat, G0 = NULL,
  ftol = 1e-8, xtol = 1e-8, gtol = 1e-8,
  max_nfev = NULL
) {
  stopifnot(
    length(lam) == 1L,
    is.finite(lam),
    lam >= 0,
    is.matrix(kernMat),
    is.numeric(kernMat),
    nrow(kernMat) >= 2L,
    nrow(kernMat) %% 2L == 0L,
    ncol(kernMat) >= 3L,
    all(is.finite(kernMat)),
    length(Gexp) == nrow(kernMat),
    all(is.finite(Gexp)),
    all(Gexp > 0),
    length(wexp) == length(Gexp),
    all(is.finite(wexp)),
    all(wexp >= 0),
    any(wexp > 0),
    length(H) == ncol(kernMat),
    all(is.finite(H))
  )

  if (!is.null(G0)) {
    stopifnot(length(G0) == 1L, is.finite(G0))
  }

  ns <- ncol(kernMat)
  start <- c(H, G0)

  fit <- .trfNoBounds(
    fun = function(x) {
      .gstarResidualLM(x, lam, Gexp, wexp, kernMat)
    },
    jac = function(x) {
      .gstarJacobianLM(x, lam, Gexp, wexp, kernMat)
    },
    x0 = start,
    ftol = ftol,
    xtol = xtol,
    gtol = gtol,
    max_nfev = max_nfev
  )

  if (!isTRUE(fit$success)) {
    stop(
      sprintf(
        paste0(
          "Frequency TRF did not converge at lambda = %.6e ",
          "(status = %d, nfev = %d)."
        ),
        lam, fit$status, fit$nfev
      ),
      call. = FALSE
    )
  }

  Hopt <- fit$x[seq_len(ns)]
  G0opt <- if (is.null(G0)) NULL else fit$x[ns + 1L]

  list(
    H = Hopt,
    G0 = G0opt,
    Gfit = .gstarKernelPrestore(Hopt, kernMat, G0opt),
    diagnostics = fit,
    solver = "experimental"
  )
}

.gstarLcurveExperimental <- function(
  Gexp, wexp, Hgs, kernMat, par, G0 = NULL
) {
  ns <- ncol(kernMat)
  npoints <- as.integer(
    par$lamDensity * (log10(par$lam_max) - log10(par$lam_min))
  )
  stopifnot(npoints >= 2L)

  ratio <- (par$lam_max / par$lam_min)^(1 / (npoints - 1L))
  lam <- par$lam_min * ratio^(0:(npoints - 1L))

  rho <- eta <- logP <- numeric(npoints)
  Hlambda <- matrix(0, ns, npoints)
  diagnostics <- vector("list", npoints)

  A <- .gstarGetAmatrix(ns)
  H <- Hgs
  logPmax <- -Inf
  first <- 1L

  for (i in rev(seq_len(npoints))) {
    fitted <- .gstarGetHExperimental(
      lam[i], Gexp, wexp, H, kernMat, G0
    )

    H <- fitted$H
    G0 <- fitted$G0
    diagnostics[i] <- list(fitted$diagnostics)

    r <- wexp * (
      1 - .gstarKernelPrestore(H, kernMat, G0) / Gexp
    )
    rho[i] <- sqrt(sum(r^2))
    eta[i] <- sqrt(sum(diff(H, differences = 2)^2))
    Hlambda[, i] <- H

    B <- .gstarGetBmatrix(H, kernMat, Gexp, wexp, G0)
    logdetC <- as.numeric(
      determinant(lam[i] * A + B, logarithm = TRUE)$modulus
    )

    # The lambda-independent prior log determinant is omitted,
    # as in the experimental time workflow.
    logP[i] <- -(rho[i]^2 + lam[i] * eta[i]^2) +
      0.5 * (ns * log(lam[i]) - logdetC) - lam[i]

    if (!is.finite(logP[i])) {
      stop(
        sprintf(
          "Non-finite frequency log evidence at lambda = %.6e.",
          lam[i]
        ),
        call. = FALSE
      )
    }

    if (logP[i] > logPmax) {
      logPmax <- logP[i]
    } else if (logP[i] < logPmax - 18) {
      first <- i
      break
    }
  }

  keep <- first:npoints
  lam <- lam[keep]
  rho <- rho[keep]
  eta <- eta[keep]
  logP <- logP[keep]
  Hlambda <- Hlambda[, keep, drop = FALSE]
  diagnostics <- diagnostics[keep]

  logP <- logP - max(logP)
  probability <- exp(logP)
  probability <- probability / sum(probability)

  lamM <- exp(sum(probability * log(lam)))

  if (par$SmFacLam > 0) {
    lamM <- exp(
      log(lamM) + par$SmFacLam * (max(log(lam)) - log(lamM))
    )
  } else if (par$SmFacLam < 0) {
    lamM <- exp(
      log(lamM) + par$SmFacLam * (log(lamM) - min(log(lam)))
    )
  }

  list(
    lamC = lamM,
    lam = lam,
    rho = rho,
    eta = eta,
    logP = logP,
    Hlambda = Hlambda,
    dH = .respectErrorBand(Hlambda, probability),
    diagnostics = diagnostics
  )
}

.gstarContSpecExperimental <- function(
  par, writeOutput = FALSE, outputDir = "output"
) {
  inp <- .gstarGetExpData(par$GstFile)
  w <- inp$w
  Gexp <- inp$Gst
  wexp <- inp$wexp
  ns <- par$ns

  wmin <- w[1L]
  wmax <- w[length(w)]

  if (par$FreqEnd == 1L) {
    smin <- exp(-pi / 2) / wmax
    smax <- exp(pi / 2) / wmin
  } else if (par$FreqEnd == 2L) {
    smin <- 1 / wmax
    smax <- 1 / wmin
  } else {
    smin <- exp(pi / 2) / wmax
    smax <- exp(-pi / 2) / wmin
  }

  ratio <- (smax / smin)^(1 / (ns - 1L))
  s <- smin * ratio^(0:(ns - 1L))
  K <- .gstarGetKernMat(s, w)

  Hraw <- -5 * rep(1, ns) + sin(pi * s)
  initial <- .gstarGetHExperimental(
    1, Gexp, wexp, Hraw, K,
    G0 = if (isTRUE(par$plateau)) min(Gexp) else NULL
  )

  scan <- NULL
  if (par$lamC == 0) {
    scan <- .gstarLcurveExperimental(
      Gexp, wexp, initial$H, K, par, initial$G0
    )
    lamC <- scan$lamC
  } else {
    lamC <- par$lamC
  }

  # Python 2.1 starts the final fit from the initialization,
  # not from the last fit in the lambda scan.
  final <- .gstarGetHExperimental(
    lamC, Gexp, wexp, initial$H, K, initial$G0
  )

  out <- list(
    s = s,
    H = final$H,
    G0 = if (is.null(final$G0)) 0.0 else final$G0,
    lamC = lamC,
    lam_C = lamC,
    w = w,
    Gexp = Gexp,
    wexp = wexp,
    kernMat = K,
    G_fit = final$Gfit,
    solver = "experimental",
    lam = if (is.null(scan)) NULL else scan$lam,
    rho = if (is.null(scan)) NULL else scan$rho,
    eta = if (is.null(scan)) NULL else scan$eta,
    log_P = if (is.null(scan)) NULL else scan$logP,
    H_lam = if (is.null(scan)) NULL else scan$Hlambda,
    dH = if (is.null(scan)) NULL else scan$dH,
    diagnostics = list(
      initial = initial$diagnostics,
      scan = if (is.null(scan)) NULL else scan$diagnostics,
      final = final$diagnostics
    ),
    reference = list(
      repository = "shane5ul/pyReSpect-freq",
      version = "2.1.0",
      commit = "908559720322941622ab39ced82e376776a39183"
    )
  )

  if (isTRUE(writeOutput)) {
    dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
    h_file <- file.path(outputDir, "H.dat")

    if (isTRUE(par$plateau)) {
      write(sprintf("# G0 = %.17e", out$G0), h_file)
      utils::write.table(
        cbind(s, out$H), h_file,
        append = TRUE, row.names = FALSE, col.names = FALSE
      )
    } else {
      utils::write.table(
        cbind(s, out$H), h_file,
        row.names = FALSE, col.names = FALSE
      )
    }

    n <- length(w)
    utils::write.table(
      cbind(w, out$G_fit[seq_len(n)], out$G_fit[n + seq_len(n)]),
      file.path(outputDir, "Gfit.dat"),
      row.names = FALSE, col.names = FALSE
    )

    if (!is.null(scan)) {
      utils::write.table(
        cbind(scan$lam, scan$rho, scan$eta),
        file.path(outputDir, "rho-eta.dat"),
        row.names = FALSE, col.names = FALSE
      )
      utils::write.table(
        cbind(scan$lam, scan$logP),
        file.path(outputDir, "logPlam.dat"),
        row.names = FALSE, col.names = FALSE
      )
    }
  }

  out
}
