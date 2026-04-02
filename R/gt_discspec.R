# gt_discspec.R — discrete relaxation spectrum, time domain (internal)
# Ported from Gt/discSpec.R; no library() / source() calls.
# Uses .gridDensity from gridDensity_shared.R and .gtGetWeights / .gtMaxwellModes.

.gtFineTuneSolution <- function(tau, t, Gexp, wexp) {
  res_tG <- function(log_tau) {
    tau_v <- exp(log_tau)
    res   <- .gtNnLLS(t, tau_v, Gexp, wexp)
    g     <- res$g
    GtM <- as.vector(exp(-outer(t, 1/tau_v)) %*% g)
    wexp * (GtM / Gexp - 1)
  }

  initErr <- sqrt(sum(res_tG(log(tau))^2))
  g       <- .gtMaxwellModes(log(tau), t, Gexp, wexp)$g
  success <- FALSE

  if (requireNamespace("minpack.lm", quietly = TRUE)) {
    tryCatch({
      fit <- minpack.lm::nls.lm(
        par     = log(tau),
        fn      = res_tG,
        control = minpack.lm::nls.lm.control(maxiter = 200)
      )
      tau_new  <- exp(fit$par)
      res_new  <- .gtMaxwellModes(log(tau_new), t, Gexp, wexp)
      finalErr <- sqrt(sum(res_tG(log(res_new$tau))^2))
      if (finalErr <= initErr) {
        tau <- res_new$tau; g <- res_new$g; success <- TRUE
      }
    }, error = function(e) invisible(NULL))
  }
  list(success = success, g = g, tau = tau)
}

.gtMergeModes <- function(g, tau, imode) {
  costFcn <- function(par) {
    gn <- par[1]; taun <- par[2]
    if (taun <= 0 || gn <= 0) return(1e10)
    g1   <- g[imode];   g2   <- g[imode + 1]
    tau1 <- tau[imode]; tau2 <- tau[imode + 1]
    tmin <- min(tau1, tau2) / 10
    tmax <- max(tau1, tau2) * 10
    integrand <- function(tt) {
      Gn <- gn * exp(-tt / taun)
      Go <- g1 * exp(-tt / tau1) + g2 * exp(-tt / tau2)
      (Gn / Go - 1)^2
    }
    pracma::integral(integrand, tmin, tmax)
  }
  ini <- c(g[imode] + g[imode + 1], 0.5 * (tau[imode] + tau[imode + 1]))
  res <- stats::optim(ini, costFcn, method = "Nelder-Mead")
  tau_new <- tau[-(imode + 1)]; tau_new[imode] <- abs(res$par[2])
  g_new   <- g[-(imode + 1)];   g_new[imode]   <- abs(res$par[1])
  list(g = g_new, tau = tau_new)
}

.gtDiscSpec <- function(par, crs, writeOutput = FALSE, outputDir = "output") {
  t    <- crs$t;  Gexp <- crs$Gexp;  wexp <- crs$wexp
  s    <- crs$s;  H    <- crs$H;     kernMat <- crs$kernMat
  n    <- length(t);  ns <- length(s)

  Nmin <- max(floor(0.5 * log10(max(t) / min(t))), 2)
  Nmax <- min(floor(3.0 * log10(max(t) / min(t))), n / 4)
  if (par$MaxNumModes > 0) Nmax <- min(Nmax, par$MaxNumModes)
  Nv   <- as.integer(seq(Nmin, Nmax))

  Gc     <- .gtKernelPrestore(H, kernMat)
  Cerror <- 1 / stats::sd(wexp * (Gc / Gexp - 1))

  wtBase <- par$deltaBaseWeightDist * seq(1, floor(1 / par$deltaBaseWeightDist) - 1)
  AICbst <- numeric(length(wtBase)); Nbst <- numeric(length(wtBase))
  nzNbst <- numeric(length(wtBase))

  loop_pb <- NULL
  if (isTRUE(par$verbose) && requireNamespace("progress", quietly = TRUE)) {
    loop_pb <- progress::progress_bar$new(
      total      = length(wtBase) * length(Nv),
      clear      = FALSE,
      show_after = 0,
      format     = "[discSpec-time] :bar :percent :current/:total"
    )
  }

  for (ib in seq_along(wtBase)) {
    wb <- wtBase[ib]
    wt <- .gtGetWeights(H, t, s, wb)
    ev <- numeric(length(Nv)); nzNv <- numeric(length(Nv))

    for (i in seq_along(Nv)) {
      if (!is.null(loop_pb)) loop_pb$tick()
      gd  <- .gridDensity(log(s), wt, Nv[i])
      res <- .gtMaxwellModes(gd$z, t, Gexp, wexp)
      ev[i]   <- res$error
      nzNv[i] <- length(res$g)
    }
    AIC        <- 2 * Nv + 2 * Cerror * ev
    AICbst[ib] <- min(AIC)
    Nbst[ib]   <- Nv[which.min(AIC)]
    nzNbst[ib] <- nzNv[which.min(AIC)]
  }

  Nopt  <- as.integer(Nbst[which.min(AICbst)])
  wbopt <- wtBase[which.min(AICbst)]

  wt  <- .gtGetWeights(H, t, s, wbopt)
  gd  <- .gridDensity(log(s), wt, Nopt)
  res <- .gtMaxwellModes(gd$z, t, Gexp, wexp)
  g   <- res$g; tau <- res$tau; error <- res$error; cKp <- res$condKp

  ft <- .gtFineTuneSolution(tau, t, Gexp, wexp)
  if (ft$success) { g <- ft$g; tau <- ft$tau }

  idx <- order(tau); tau <- tau[idx]; g <- g[idx]
  if (length(tau) > 1) {
    tauSpacing <- tau[-1] / tau[-length(tau)]
    itry <- 0
    while (min(tauSpacing) < par$minTauSpacing && itry < 3) {
      imode <- which.min(tauSpacing)
      mg    <- .gtMergeModes(g, tau, imode)
      g <- mg$g; tau <- mg$tau
      ft <- .gtFineTuneSolution(tau, t, Gexp, wexp)
      if (ft$success) { g <- ft$g; tau <- ft$tau }
      if (length(tau) > 1) tauSpacing <- tau[-1] / tau[-length(tau)] else break
      itry <- itry + 1
    }
  }

  if (isTRUE(writeOutput)) {
    if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)
    utils::write.table(cbind(g, tau),               file.path(outputDir, "dmodes.dat"),
                       row.names = FALSE, col.names = FALSE)
    utils::write.table(cbind(wtBase, nzNbst, AICbst), file.path(outputDir, "aic.dat"),
                       row.names = FALSE, col.names = FALSE)
  }

  if (par$verbose) {
    cat(sprintf("(*) Nopt = %d\n", length(g)))
    cat(sprintf("(*) log10(Condition number) = %.2f\n", log10(cKp)))
    cat("\n  i\t\tg(i)\t\ttau(i)\n-------------------------------------------\n")
    for (i in seq_along(g))
      cat(sprintf("%3d\t%.5e\t%.5e\n", i, g[i], tau[i]))
  }

  if (par$plotting && isTRUE(writeOutput)) {
    grDevices::pdf(file.path(outputDir, "dmodes_time.pdf"))
    graphics::plot(tau, g, log = "xy", type = "o", col = "blue",
                   xlab = expression(tau), ylab = "g",
                   main = "Discrete Relaxation Spectrum (time)")
    graphics::lines(s, exp(H), col = "red", lwd = 2)
    graphics::legend("bottomright", c("discrete", "continuous"),
                     col = c("blue", "red"), lty = 1)
    grDevices::dev.off()
  }

  list(g = g, tau = tau, error = error, condKp = cKp, Nopt = Nopt)
}
