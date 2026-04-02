# gstar_discspec.R — discrete relaxation spectrum, frequency domain (internal)
# Ported from Gstar/discSpec.R; no library() / source() calls.
# Uses .gridDensity from gridDensity_shared.R and .gstar* from gstar_common.R.

# ── nnLLS (frequency domain, optional plateau column) ─────────────────────────
.gstarNnLLS <- function(w, tau, Gexp, wexp, isPlateau) {
  n    <- length(Gexp) / 2
  ntau <- length(tau)

  ws  <- outer(w, tau)           # n×ntau: ws[i,j] = w[i]*tau[j]
  ws2 <- ws^2
  K   <- rbind(ws2 / (1 + ws2), ws / (1 + ws2))    # 2n x ntau

  if (isPlateau) {
    col_G0 <- c(rep(1, n), rep(0, n))
    K      <- cbind(K, col_G0)
  }

  Kp     <- (wexp / Gexp) * K
  condKp <- kappa(Kp, exact = FALSE)
  g      <- nnls::nnls(Kp, wexp)$x

  GstM  <- as.vector(K %*% g)
  error <- sum((wexp * (GstM / Gexp - 1))^2)
  list(g = g, error = error, condKp = condKp)
}

# ── MaxwellModes (frequency domain, with plateau support) ─────────────────────
.gstarDiscMaxwellModes <- function(z, w, Gexp, wexp, isPlateau) {
  z   <- pmax(z, log(0.02 / max(w)))
  z   <- pmin(z, log(50  / min(w)))
  tau <- exp(z)
  n   <- length(w)

  res    <- .gstarNnLLS(w, tau, Gexp, wexp, isPlateau)
  g      <- res$g;  error <- res$error;  condKp <- res$condKp

  izero <- which(max(w) * tau < 0.02 | min(w) * tau > 50)
  if (length(izero) > 0) { tau <- tau[-izero]; g <- g[-izero] }

  if (isPlateau) {
    g_main <- g[seq_len(length(g) - 1)]
    if (length(g_main) > 0 && max(g_main) > 0)
      izero <- which(g_main / max(g_main) < 1e-8)
    else
      izero <- integer(0)
  } else {
    if (length(g) > 0 && max(g) > 0)
      izero <- which(g / max(g) < 1e-8)
    else
      izero <- integer(0)
  }
  if (length(izero) > 0) { tau <- tau[-izero]; g <- g[-izero] }

  list(g = g, tau = tau, error = error, condKp = condKp)
}

# ── GetWeights (frequency domain) ─────────────────────────────────────────────
.gstarGetWeights <- function(H, w, s, wb) {
  ns <- length(s);  n <- length(w)

  hs <- numeric(ns)
  hs[1]        <- 0.5 * log(s[2] / s[1])
  hs[ns]       <- 0.5 * log(s[ns] / s[ns - 1])
  hs[2:(ns-1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns - 2)]))

  ws   <- outer(w, s)                               # n×ns
  ws2  <- ws^2
  kern <- rbind(ws2 / (1 + ws2), ws / (1 + ws2))   # 2n x ns
  scale <- hs * exp(H)
  K     <- as.vector(kern %*% scale)
  wij   <- sweep(kern, 2, scale, "*")               # column-scale
  wij   <- sweep(wij,  1, K,    "/")               # row-normalise

  wt <- colSums(wij)
  wt <- wt / pracma::trapz(log(s), wt)
  wt <- (1 - wb) * wt + wb * mean(wt) * rep(1, ns)
  wt
}

# ── mergeModes_magic ───────────────────────────────────────────────────────────
.gstarMergeModes <- function(g, tau, imode) {
  i1 <- imode;  i2 <- imode + 1
  g1 <- g[i1];  g2 <- g[i2];  tau1 <- tau[i1];  tau2 <- tau[i2]
  wmin <- min(1 / tau1, 1 / tau2) / 10
  wmax <- max(1 / tau1, 1 / tau2) * 10

  normKern <- function(ww, gn, taun) {
    wt  <- ww * taun; Gnp  <- gn * (wt^2) / (1 + wt^2); Gnpp <- gn * wt / (1 + wt^2)
    wt  <- ww * tau1; Gop  <- g1 * (wt^2) / (1 + wt^2); Gopp <- g1 * wt / (1 + wt^2)
    wt  <- ww * tau2; Gop  <- Gop  + g2 * (wt^2) / (1 + wt^2)
                      Gopp <- Gopp + g2 * wt / (1 + wt^2)
    (Gnp / Gop - 1)^2 + (Gnpp / Gopp - 1)^2
  }

  costFcn <- function(par) {
    stats::integrate(normKern, wmin, wmax, gn = par[1], taun = par[2])$value
  }

  ini <- c(g1 + g2, 0.5 * (tau1 + tau2))
  res <- stats::optim(ini, costFcn, method = "Nelder-Mead")

  newtau        <- tau[-i2];  newtau[i1] <- res$par[2]
  newg          <- g[-i2];    newg[i1]   <- res$par[1]
  list(g = newg, tau = newtau)
}

# ── FineTuneSolution ───────────────────────────────────────────────────────────
.gstarFineTuneSolution <- function(tau, w, Gexp, wexp, isPlateau) {
  res_wG <- function(tau_vec) {
    g   <- .gstarNnLLS(w, tau_vec, Gexp, wexp, isPlateau)$g
    n   <- length(w)
    ws  <- outer(w, tau_vec)
    ws2 <- ws^2
    K   <- rbind(ws2 / (1 + ws2), ws / (1 + ws2))
    if (isPlateau) {
      Gmodel         <- as.vector(K %*% g[-length(g)])
      Gmodel[1:n]    <- Gmodel[1:n] + g[length(g)]
    } else {
      Gmodel <- as.vector(K %*% g)
    }
    wexp * (Gmodel / Gexp - 1)
  }

  success   <- FALSE
  initError <- sqrt(sum(res_wG(tau)^2))

  tau_new <- tryCatch({
    res <- minpack.lm::nls.lm(
      par     = tau,
      fn      = res_wG,
      lower   = rep(0.02 / max(w), length(tau)),
      upper   = rep(50  / min(w),  length(tau)),
      control = minpack.lm::nls.lm.control(maxiter = 200)
    )
    res$par
  }, error = function(e) tau)

  res_MM  <- .gstarDiscMaxwellModes(log(tau_new), w, Gexp, wexp, isPlateau)
  tau_out <- res_MM$tau

  finalError <- sqrt(sum(res_wG(tau_out)^2))
  if (finalError < initError) success <- TRUE else tau_out <- tau

  res_MM2 <- .gstarDiscMaxwellModes(log(tau_out), w, Gexp, wexp, isPlateau)
  list(success = success, g = res_MM2$g, tau = res_MM2$tau)
}

# ── initializeDiscSpec ─────────────────────────────────────────────────────────
.gstarInitializeDiscSpec <- function(par, contResult) {
  if (is.null(contResult) || is.null(contResult$s) || is.null(contResult$H))
    stop("contResult must include s and H for frequency-domain discrete solve")

  if (par$verbose)
    cat(sprintf("\n(*) Start\n(*) Loading Data Files: %s...\n", par$GstFile))

  inp  <- .gstarGetExpData(par$GstFile)
  w    <- inp$w;  Gexp <- inp$Gst;  wexp <- inp$wexp
  n    <- length(w)

  s  <- contResult$s
  H  <- contResult$H
  G0 <- if (!is.null(contResult$G0)) contResult$G0 else 0

  ns   <- length(s)
  Nmax <- min(floor(3.0 * log10(max(w) / min(w))), n / 4)
  if (par$MaxNumModes > 0) Nmax <- min(Nmax, par$MaxNumModes)
  Nmin <- max(floor(0.5 * log10(max(w) / min(w))), 3)
  Nv   <- as.integer(seq(Nmin, Nmax))

  kernMat <- .gstarGetKernMat(s, w)
  Gc      <- if (isTRUE(par$plateau)) {
    .gstarKernelPrestore(H, kernMat, G0)
  } else {
    .gstarKernelPrestore(H, kernMat)
  }
  Cerror <- 1 / stats::sd(wexp * (Gc / Gexp - 1))

  list(w = w, Gexp = Gexp, wexp = wexp, s = s, H = H, G0 = G0,
       Nv = Nv, Gc = Gc, Cerror = Cerror)
}

# ── getDiscSpecMagic (main entry point) ───────────────────────────────────────
.gstarGetDiscSpecMagic <- function(par, writeOutput = FALSE, outputDir = "output",
                                   contResult = NULL) {
  res0   <- .gstarInitializeDiscSpec(par, contResult)
  w      <- res0$w;   Gexp <- res0$Gexp;  wexp  <- res0$wexp
  s      <- res0$s;   H    <- res0$H;     Nv    <- res0$Nv
  Gc     <- res0$Gc;  Cerror <- res0$Cerror
  n      <- length(w)

  wtBase <- par$deltaBaseWeightDist * seq(1, floor(1 / par$deltaBaseWeightDist) - 1)
  AICbst <- numeric(length(wtBase))
  Nbst   <- numeric(length(wtBase))
  nzNbst <- numeric(length(wtBase))
  npts   <- length(Nv)

  loop_pb <- NULL
  if (isTRUE(par$verbose) && requireNamespace("progress", quietly = TRUE)) {
    loop_pb <- progress::progress_bar$new(
      total      = length(wtBase) * npts,
      clear      = FALSE,
      show_after = 0,
      format     = "[discSpec-freq] :bar :percent :current/:total"
    )
  }

  for (ib in seq_along(wtBase)) {
    wb   <- wtBase[ib]
    wt   <- .gstarGetWeights(H, w, s, wb)
    ev   <- numeric(npts)
    nzNv <- numeric(npts)

    for (i in seq_along(Nv)) {
      if (!is.null(loop_pb)) loop_pb$tick()
      gd      <- .gridDensity(log(s), wt, Nv[i])
      mm      <- .gstarDiscMaxwellModes(gd$z, w, Gexp, wexp, isTRUE(par$plateau))
      ev[i]   <- mm$error
      nzNv[i] <- length(mm$g)
    }

    AIC        <- 2 * Nv + 2 * Cerror * ev
    AICbst[ib] <- min(AIC)
    Nbst[ib]   <- Nv[which.min(AIC)]
    nzNbst[ib] <- nzNv[which.min(AIC)]
  }

  Nopt  <- as.integer(Nbst[which.min(AICbst)])
  wbopt <- wtBase[which.min(AICbst)]

  wt  <- .gstarGetWeights(H, w, s, wbopt)
  gd  <- .gridDensity(log(s), wt, Nopt)
  mm  <- .gstarDiscMaxwellModes(gd$z, w, Gexp, wexp, isTRUE(par$plateau))
  g   <- mm$g;  tau <- mm$tau;  error <- mm$error;  cKp <- mm$condKp

  ft <- .gstarFineTuneSolution(tau, w, Gexp, wexp, isTRUE(par$plateau))
  if (ft$success) { g <- ft$g;  tau <- ft$tau }

  idx        <- order(tau)
  tau        <- tau[idx]
  if (isTRUE(par$plateau)) {
    g[seq_len(length(g) - 1)] <- g[idx]
  } else {
    g <- g[idx]
  }

  tauSpacing <- if (length(tau) > 1) tau[-1] / tau[-length(tau)] else Inf
  itry <- 0
  while (length(tauSpacing) > 0 && min(tauSpacing) < par$minTauSpacing && itry < 3) {
    cat("\tTau Spacing < minTauSpacing\n")
    imode      <- which.min(tauSpacing)
    mg         <- .gstarMergeModes(g, tau, imode)
    g <- mg$g;  tau <- mg$tau
    ft <- .gstarFineTuneSolution(tau, w, Gexp, wexp, isTRUE(par$plateau))
    if (ft$success) { g <- ft$g;  tau <- ft$tau }
    tauSpacing <- if (length(tau) > 1) tau[-1] / tau[-length(tau)] else Inf
    itry <- itry + 1
  }

  G0 <- 0
  if (isTRUE(par$plateau)) { G0 <- g[length(g)];  g <- g[-length(g)] }

  if (isTRUE(writeOutput)) {
    if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)
  }

  if (par$verbose) {
    cat(sprintf("(*) Number of optimum nodes = %d\n", length(g)))
    cat(sprintf("(*) Condition number: %e\n", cKp))
    if (isTRUE(par$plateau)) cat(sprintf("(*) Plateau Modulus G0 = %.3e\n", G0))
  }

  if (isTRUE(writeOutput)) {
    if (isTRUE(par$plateau)) {
      write(sprintf("# G0 = %e", G0), file.path(outputDir, "dmodes.dat"))
      utils::write.table(cbind(g, tau), file.path(outputDir, "dmodes.dat"),
                         append = TRUE, row.names = FALSE, col.names = FALSE)
    } else {
      utils::write.table(cbind(g, tau), file.path(outputDir, "dmodes.dat"),
                         row.names = FALSE, col.names = FALSE)
    }
    utils::write.table(cbind(wtBase, nzNbst, AICbst), file.path(outputDir, "aic.dat"),
                       row.names = FALSE, col.names = FALSE)

  ws  <- outer(w, tau)
  ws2 <- ws^2
  K   <- rbind(ws2 / (1 + ws2), ws / (1 + ws2))
    GstM <- as.vector(K %*% g)
    if (isTRUE(par$plateau)) GstM[1:n] <- GstM[1:n] + G0
    utils::write.table(cbind(w, GstM[1:n], GstM[(n+1):(2*n)]),
                       file.path(outputDir, "Gfitd.dat"),
                       row.names = FALSE, col.names = FALSE)
  }

  if (par$verbose) {
    cat("\n\t\tModes\n\t\t-----\n\n")
    cat("  i \t    g(i) \t    tau(i)\n")
    cat("-------------------------------------------\n")
    for (i in seq_along(g))
      cat(sprintf("%3d \t %.5e \t %.5e\n", i, g[i], tau[i]))
  }

  list(Nopt = Nopt, g = g, tau = tau, error = error)
}
