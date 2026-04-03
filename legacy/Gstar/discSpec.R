# discSpec.R
# Traduzione fedele di discSpec.py (pyReSpect-freq)

library(minpack.lm)
library(pracma)  # per interp1, cumtrapz

# ── nnLLS ─────────────────────────────────────────────────────────────────────
nnLLS <- function(w, tau, Gexp, wexp, isPlateau) {
  n    <- length(Gexp) / 2
  ntau <- length(tau)
  res  <- meshgrid(tau, w)
  S <- res$X; W <- res$Y
  ws  <- S * W; ws2 <- ws^2
  K   <- rbind(ws2/(1+ws2), ws/(1+ws2))  # 2n x ntau

  if (isPlateau) {
    col_G0      <- c(rep(1, n), rep(0, n))
    K           <- cbind(K, col_G0)
  }

  Kp     <- (wexp / Gexp) * K   # broadcasting row-wise (equiv. np.diag(wexp/Gexp) @ K)
  condKp <- kappa(Kp, exact = FALSE)
  g      <- nnls::nnls(Kp, wexp)$x

  GstM  <- as.vector(K %*% g)
  error <- sum((wexp * (GstM/Gexp - 1))^2)

  return(list(g = g, error = error, condKp = condKp))
}

# ── MaxwellModes ──────────────────────────────────────────────────────────────
MaxwellModes <- function(z, w, Gexp, wexp, isPlateau) {
  # clamp z per evitare tau fuori dalla finestra fisica
  z   <- pmax(z, log(0.02/max(w)))
  z   <- pmin(z, log(50/min(w)))
  tau <- exp(z)
  n   <- length(w)

  res <- nnLLS(w, tau, Gexp, wexp, isPlateau)
  g   <- res$g; error <- res$error; condKp <- res$condKp

  # rimuovi modi fuori finestra con peso grande
  izero <- which(max(w)*tau < 0.02 | min(w)*tau > 50)
  if (length(izero) > 0) { tau <- tau[-izero]; g <- g[-izero] }

  # rimuovi pesi piccoli
  if (isPlateau) {
    g_main <- g[seq_len(length(g)-1)]
    izero  <- which(g_main / max(g_main) < 1e-8)
  } else {
    izero  <- which(g / max(g) < 1e-8)
  }
  if (length(izero) > 0) { tau <- tau[-izero]; g <- g[-izero] }

  return(list(g = g, tau = tau, error = error, condKp = condKp))
}

# ── GetWeights ────────────────────────────────────────────────────────────────
GetWeights <- function(H, w, s, wb) {
  ns <- length(s); n <- length(w)
  hs <- numeric(ns)
  hs[1]        <- 0.5 * log(s[2]/s[1])
  hs[ns]       <- 0.5 * log(s[ns]/s[ns-1])
  hs[2:(ns-1)] <- 0.5 * (log(s[3:ns]) - log(s[1:(ns-2)]))

  res <- meshgrid(s, w)
  S <- res$X; W <- res$Y
  ws <- S*W; ws2 <- ws^2
  kern <- rbind(ws2/(1+ws2), ws/(1+ws2))  # 2n x ns

  wij <- kern * matrix(rep(hs * exp(H), each=2*n), nrow=2*n)  # 2n x ns
  K   <- as.vector(kern %*% (hs * exp(H)))                    # 2n x 1

  for (i in seq_len(n)) wij[i,] <- wij[i,] / K[i]

  wt <- numeric(ns)
  for (j in seq_len(ns)) wt[j] <- sum(wij[,j])

  wt <- wt / trapz(log(s), wt)
  wt <- (1 - wb) * wt + (wb * mean(wt)) * rep(1, ns)
  return(wt)
}

# ── GridDensity ───────────────────────────────────────────────────────────────
GridDensity <- function(x, px, N) {
  npts <- 100
  xi   <- seq(min(x), max(x), length.out = npts)
  pint <- interp1(x, px, xi, method = "spline")
  ci   <- cumtrapz(xi, pint)
  pint <- pint / ci[npts]
  ci   <- ci   / ci[npts]

  alfa <- 1 / (N - 1)
  z    <- numeric(N)
  z[1] <- min(x); z[N] <- max(x)

  beta     <- seq(0.5, N-0.5) * alfa
  fint_inv <- approxfun(ci, xi, rule = 2)
  zij      <- c(z[1], fint_inv(beta[seq_len(N-1)]), z[N])
  h        <- diff(zij)

  beta2    <- seq(1, N-1) * alfa
  z[2:(N-1)] <- fint_inv(beta2[seq_len(N-2)])

  return(list(z = z, h = h))
}

# ── mergeModes_magic ──────────────────────────────────────────────────────────
mergeModes_magic <- function(g, tau, imode) {
  i1 <- imode; i2 <- imode + 1
  g1 <- g[i1]; g2 <- g[i2]; tau1 <- tau[i1]; tau2 <- tau[i2]
  wmin <- min(1/tau1, 1/tau2) / 10
  wmax <- max(1/tau1, 1/tau2) * 10

  normKern <- function(w, gn, taun) {
    wt <- w*taun; Gnp  <- gn*(wt^2)/(1+wt^2); Gnpp <- gn*wt/(1+wt^2)
    wt <- w*tau1; Gop  <- g1*(wt^2)/(1+wt^2); Gopp <- g1*wt/(1+wt^2)
    wt <- w*tau2; Gop  <- Gop  + g2*(wt^2)/(1+wt^2)
                  Gopp <- Gopp + g2*wt/(1+wt^2)
    (Gnp/Gop - 1)^2 + (Gnpp/Gopp - 1)^2
  }

  costFcn <- function(par) {
    integrate(normKern, wmin, wmax, gn=par[1], taun=par[2])$value
  }

  ini  <- c(g1 + g2, 0.5*(tau1+tau2))
  res  <- optim(ini, costFcn, method = "Nelder-Mead")

  newtau        <- tau[-i2]; newtau[i1] <- res$par[2]
  newg          <- g[-i2];   newg[i1]   <- res$par[1]
  return(list(g = newg, tau = newtau))
}

# ── FineTuneSolution ──────────────────────────────────────────────────────────
FineTuneSolution <- function(tau, w, Gexp, wexp, isPlateau) {
  res_wG <- function(tau_vec) {
    g <- nnLLS(w, tau_vec, Gexp, wexp, isPlateau)$g
    n <- length(w)
    res2 <- meshgrid(tau_vec, w)
    S <- res2$X; W <- res2$Y
    ws <- S*W; ws2 <- ws^2
    K  <- rbind(ws2/(1+ws2), ws/(1+ws2))
    if (isPlateau) {
      Gmodel <- as.vector(K %*% g[-length(g)])
      Gmodel[1:n] <- Gmodel[1:n] + g[length(g)]
    } else {
      Gmodel <- as.vector(K %*% g)
    }
    wexp * (Gmodel/Gexp - 1)
  }

  success   <- FALSE
  initError <- sqrt(sum(res_wG(tau)^2))
  n         <- length(w)

  tau_new <- tryCatch({
    res <- nls.lm(par = tau,
                  fn  = res_wG,
                  lower = rep(0.02/max(w), length(tau)),
                  upper = rep(50/min(w),   length(tau)),
                  control = nls.lm.control(maxiter=200))
    res$par
  }, error = function(e) tau)

  res_MM <- MaxwellModes(log(tau_new), w, Gexp, wexp, isPlateau)
  g <- res_MM$g; tau_out <- res_MM$tau

  finalError <- sqrt(sum(res_wG(tau_out)^2))
  if (finalError < initError) success <- TRUE else tau_out <- tau

  res_MM2 <- MaxwellModes(log(tau_out), w, Gexp, wexp, isPlateau)
  return(list(success = success, g = res_MM2$g, tau = res_MM2$tau))
}

# ── initializeDiscSpec ────────────────────────────────────────────────────────
initializeDiscSpec <- function(par) {
  if (par$verbose)
    cat(sprintf("\n(*) Start\n(*) Loading Data Files: %s...\n", par$GstFile))

  inp  <- GetExpData(par$GstFile)
  w    <- inp$w; Gexp <- inp$Gst; wexp <- inp$wexp
  n    <- length(w)

  # Leggi H.dat
  hlines <- readLines("output/H.dat")
  hlines <- hlines[!grepl("^#", hlines) & nzchar(trimws(hlines))]
  hdat   <- read.table(text = paste(hlines, collapse="\n"), header=FALSE)
  s <- hdat[[1]]; H <- hdat[[2]]
  ns <- length(s)

  Nmax <- min(floor(3.0 * log10(max(w)/min(w))), n/4)
  if (par$MaxNumModes > 0) Nmax <- min(Nmax, par$MaxNumModes)
  Nmin <- max(floor(0.5 * log10(max(w)/min(w))), 3)
  Nv   <- as.integer(seq(Nmin, Nmax))

  kernMat <- getKernMat(s, w)

  G0 <- 0
  if (isTRUE(par$plateau)) {
    hdr <- readLines("output/H.dat")[1]
    G0  <- as.numeric(regmatches(hdr, regexpr("[0-9eE+\\-\\.]+$", hdr)))
    if (is.na(G0)) G0 <- 0
    Gc  <- kernel_prestore(H, kernMat, G0)
  } else {
    Gc  <- kernel_prestore(H, kernMat)
  }

  Cerror <- 1 / sd(wexp * (Gc/Gexp - 1))

  return(list(w=w, Gexp=Gexp, wexp=wexp, s=s, H=H, Nv=Nv, Gc=Gc, Cerror=Cerror))
}

# ── getDiscSpecMagic ──────────────────────────────────────────────────────────
getDiscSpecMagic <- function(par) {
  res0 <- initializeDiscSpec(par)
  w <- res0$w; Gexp <- res0$Gexp; wexp <- res0$wexp
  s <- res0$s; H <- res0$H; Nv <- res0$Nv
  Gc <- res0$Gc; Cerror <- res0$Cerror
  n <- length(w)

  wtBase <- par$deltaBaseWeightDist * seq(1, floor(1/par$deltaBaseWeightDist)-1)
  AICbst  <- numeric(length(wtBase))
  Nbst    <- numeric(length(wtBase))
  nzNbst  <- numeric(length(wtBase))
  npts    <- length(Nv)

  for (ib in seq_along(wtBase)) {
    wb  <- wtBase[ib]
    wt  <- GetWeights(H, w, s, wb)
    ev  <- numeric(npts)
    nzNv<- numeric(npts)

    for (i in seq_along(Nv)) {
      N      <- Nv[i]
      gd     <- GridDensity(log(s), wt, N)
      mm     <- MaxwellModes(gd$z, w, Gexp, wexp, isTRUE(par$plateau))
      ev[i]  <- mm$error
      nzNv[i]<- length(mm$g)
    }

    AIC        <- 2*Nv + 2*Cerror*ev
    AICbst[ib] <- min(AIC)
    Nbst[ib]   <- Nv[which.min(AIC)]
    nzNbst[ib] <- nzNv[which.min(AIC)]
  }

  Nopt  <- as.integer(Nbst[which.min(AICbst)])
  wbopt <- wtBase[which.min(AICbst)]

  wt             <- GetWeights(H, w, s, wbopt)
  gd             <- GridDensity(log(s), wt, Nopt)
  mm             <- MaxwellModes(gd$z, w, Gexp, wexp, isTRUE(par$plateau))
  g <- mm$g; tau <- mm$tau; error <- mm$error; cKp <- mm$condKp

  ft <- FineTuneSolution(tau, w, Gexp, wexp, isTRUE(par$plateau))
  if (ft$success) { g <- ft$g; tau <- ft$tau }

  # Ordina tau e mergia modi vicini
  idx        <- order(tau)
  tau        <- tau[idx]
  if (isTRUE(par$plateau)) g[seq_len(length(g)-1)] <- g[idx] else g <- g[idx]
  tauSpacing <- tau[-1] / tau[-length(tau)]
  itry       <- 0

  while (length(tauSpacing) > 0 && min(tauSpacing) < par$minTauSpacing && itry < 3) {
    cat("\tTau Spacing < minTauSpacing\n")
    imode      <- which.min(tauSpacing)
    mg         <- mergeModes_magic(g, tau, imode)
    g <- mg$g; tau <- mg$tau
    ft <- FineTuneSolution(tau, w, Gexp, wexp, isTRUE(par$plateau))
    if (ft$success) { g <- ft$g; tau <- ft$tau }
    tauSpacing <- tau[-1] / tau[-length(tau)]
    itry <- itry + 1
  }

  G0 <- 0
  if (isTRUE(par$plateau)) { G0 <- g[length(g)]; g <- g[-length(g)] }

  if (par$verbose) {
    cat(sprintf("(*) Number of optimum nodes = %d\n", length(g)))
    cat(sprintf("(*) Condition number: %e\n", cKp))
    if (isTRUE(par$plateau)) cat(sprintf("(*) Plateau Modulus G0 = %.3e\n", G0))

    if (isTRUE(par$plateau)) {
      write(sprintf("# G0 = %e", G0), "output/dmodes.dat")
      write.table(cbind(g, tau), "output/dmodes.dat", append=TRUE,
                  row.names=FALSE, col.names=FALSE)
    } else {
      write.table(cbind(g, tau), "output/dmodes.dat", row.names=FALSE, col.names=FALSE)
    }

    write.table(cbind(wtBase, nzNbst, AICbst), "output/aic.dat",
                row.names=FALSE, col.names=FALSE)

    # Gfitd.dat
    res2 <- meshgrid(tau, w)
    S <- res2$X; W2 <- res2$Y
    ws <- S*W2; ws2 <- ws^2
    K  <- rbind(ws2/(1+ws2), ws/(1+ws2))
    GstM <- as.vector(K %*% g)
    if (isTRUE(par$plateau)) GstM[1:n] <- GstM[1:n] + G0
    write.table(cbind(w, GstM[1:n], GstM[(n+1):(2*n)]),
                "output/Gfitd.dat", row.names=FALSE, col.names=FALSE)

    cat("\n\t\tModes\n\t\t-----\n\n")
    cat("  i \t    g(i) \t    tau(i)\n")
    cat("-------------------------------------------\n")
    for (i in seq_along(g))
      cat(sprintf("%3d \t %.5e \t %.5e\n", i, g[i], tau[i]))
  }

  return(list(Nopt=Nopt, g=g, tau=tau, error=error))
}
