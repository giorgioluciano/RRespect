
# discSpec_time.R
# Faithfully translated from discSpec.py (pyReSpect-time) getDiscSpecMagic()
# mergeModes uses exp(-t/tau) kernel (not Maxwell freq kernel)

library(pracma);  library(nnls)
source("common.R")
source("GetWeights.R")
source("GridDensity.R")
source("MaxwellModes.R")

FineTuneSolution_time <- function(tau, t, Gexp, wexp) {
  res_tG <- function(log_tau) {
    tau_v <- exp(log_tau)
    res   <- nnLLS_time(t, tau_v, Gexp, wexp)
    g     <- res$g
    res2  <- meshgrid(tau_v, t)
    S <- res2$X; T2 <- res2$Y
    GtM <- as.vector(exp(-T2/S) %*% g)
    return(wexp * (GtM / Gexp - 1))
  }

  initErr <- sqrt(sum(res_tG(log(tau))^2))
  g    <- MaxwellModes(log(tau), t, Gexp, wexp)$g
  success <- FALSE

  if (requireNamespace("minpack.lm", quietly = TRUE)) {
    tryCatch({
      fit <- minpack.lm::nls.lm(
        par     = log(tau),
        fn      = res_tG,
        control = minpack.lm::nls.lm.control(maxiter = 200)
      )
      tau_new  <- exp(fit$par)
      res_new  <- MaxwellModes(log(tau_new), t, Gexp, wexp)
      finalErr <- sqrt(sum(res_tG(log(res_new$tau))^2))
      if (finalErr <= initErr) {
        tau <- res_new$tau; g <- res_new$g; success <- TRUE
      }
    }, error = function(e) invisible(NULL))
  }
  return(list(success = success, g = g, tau = tau))
}

mergeModes_time <- function(g, tau, imode) {
  costFcn <- function(par) {
    gn <- par[1]; taun <- par[2]
    if (taun <= 0 || gn <= 0) return(1e10)
    g1 <- g[imode]; g2 <- g[imode+1]
    tau1 <- tau[imode]; tau2 <- tau[imode+1]
    tmin <- min(tau1, tau2) / 10
    tmax <- max(tau1, tau2) * 10
    integrand <- function(tt) {
      Gn <- gn * exp(-tt/taun)
      Go <- g1 * exp(-tt/tau1) + g2 * exp(-tt/tau2)
      (Gn/Go - 1)^2
    }
    pracma::integral(integrand, tmin, tmax)
  }
  ini <- c(g[imode]+g[imode+1], 0.5*(tau[imode]+tau[imode+1]))
  res <- optim(ini, costFcn, method = "Nelder-Mead")
  tau_new <- tau[-(imode+1)]; tau_new[imode] <- abs(res$par[2])
  g_new   <- g[-(imode+1)];   g_new[imode]   <- abs(res$par[1])
  return(list(g = g_new, tau = tau_new))
}

discSpec <- function(par, crs) {
  t    <- crs$t;  Gexp <- crs$Gexp;  wexp <- crs$wexp
  s    <- crs$s;  H    <- crs$H;     kernMat <- crs$kernMat
  n    <- length(t);  ns <- length(s)

  Nmin <- max(floor(0.5 * log10(max(t)/min(t))), 2)
  Nmax <- min(floor(3.0 * log10(max(t)/min(t))), n/4)
  if (par$MaxNumModes > 0) Nmax <- min(Nmax, par$MaxNumModes)
  Nv   <- as.integer(seq(Nmin, Nmax))

  Gc     <- kernel_prestore(H, kernMat)
  Cerror <- 1 / sd(wexp * (Gc / Gexp - 1))

  wtBase <- par$deltaBaseWeightDist * seq(1, floor(1/par$deltaBaseWeightDist) - 1)
  AICbst <- numeric(length(wtBase)); Nbst <- numeric(length(wtBase))
  nzNbst <- numeric(length(wtBase))

  for (ib in seq_along(wtBase)) {
    wb <- wtBase[ib]
    wt <- GetWeights(H, t, s, wb)
    ev <- numeric(length(Nv)); nzNv <- numeric(length(Nv))

    for (i in seq_along(Nv)) {
      gd  <- GridDensity(log(s), wt, Nv[i])
      res <- MaxwellModes(gd$z, t, Gexp, wexp)
      ev[i]   <- res$error
      nzNv[i] <- length(res$g)
    }
    AIC        <- 2*Nv + 2*Cerror*ev
    AICbst[ib] <- min(AIC)
    Nbst[ib]   <- Nv[which.min(AIC)]
    nzNbst[ib] <- nzNv[which.min(AIC)]
  }

  Nopt  <- as.integer(Nbst[which.min(AICbst)])
  wbopt <- wtBase[which.min(AICbst)]

  wt    <- GetWeights(H, t, s, wbopt)
  gd    <- GridDensity(log(s), wt, Nopt)
  res   <- MaxwellModes(gd$z, t, Gexp, wexp)
  g     <- res$g; tau <- res$tau; error <- res$error; cKp <- res$condKp

  ft <- FineTuneSolution_time(tau, t, Gexp, wexp)
  if (ft$success) { g <- ft$g; tau <- ft$tau }

  idx <- order(tau); tau <- tau[idx]; g <- g[idx]
  if (length(tau) > 1) {
    tauSpacing <- tau[-1] / tau[-length(tau)]
    itry <- 0
    while (min(tauSpacing) < par$minTauSpacing && itry < 3) {
      imode <- which.min(tauSpacing)
      mg    <- mergeModes_time(g, tau, imode)
      g <- mg$g; tau <- mg$tau
      ft <- FineTuneSolution_time(tau, t, Gexp, wexp)
      if (ft$success) { g <- ft$g; tau <- ft$tau }
      if (length(tau) > 1) tauSpacing <- tau[-1]/tau[-length(tau)] else break
      itry <- itry + 1
    }
  }

  if (!dir.exists("output")) dir.create("output")
  write.table(cbind(g, tau), "output/dmodes.dat",
              row.names=FALSE, col.names=FALSE)
  write.table(cbind(wtBase, nzNbst, AICbst), "output/aic.dat",
              row.names=FALSE, col.names=FALSE)

  if (par$verbose) {
    cat(sprintf("(*) Nopt = %d\n", length(g)))
    cat(sprintf("(*) log10(Condition number) = %.2f\n", log10(cKp)))
    cat("\n  i\t\tg(i)\t\ttau(i)\n-------------------------------------------\n")
    for (i in seq_along(g))
      cat(sprintf("%3d\t%.5e\t%.5e\n", i, g[i], tau[i]))
  }

  if (par$plotting) {
    pdf("output/dmodes_time.pdf")
    plot(tau, g, log="xy", type="o", col="blue",
         xlab=expression(tau), ylab="g", main="Discrete Relaxation Spectrum (time)")
    lines(s, exp(H), col="red", lwd=2)
    legend("bottomright", c("discrete","continuous"), col=c("blue","red"), lty=1)
    dev.off()
  }

  return(list(g = g, tau = tau, error = error, condKp = cKp, Nopt = Nopt))
}
