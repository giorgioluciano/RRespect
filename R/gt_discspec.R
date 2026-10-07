# gt_discspec.R — discrete relaxation spectrum, time domain (internal)
# Ported from Gt/discSpec.R; no library() / source() calls.
# Uses .gridDensity from gridDensity_shared.R and .gtGetWeights / .gtMaxwellModes.

.gtFineTuneSolution <- function(
  tau, t, Gexp, wexp, isPlateau = FALSE
) {
  res_tG <- function(log_tau) {
    tau_v <- exp(log_tau)
    res <- .gtNnLLS(
      t, tau_v, Gexp, wexp, isPlateau
    )
    g <- res$g

    K <- exp(-outer(t, 1 / tau_v))

    if (isTRUE(isPlateau)) {
      GtM <- as.vector(K %*% g[seq_along(tau_v)]) +
        g[length(g)]
    } else {
      GtM <- as.vector(K %*% g)
    }

    wexp * (GtM / Gexp - 1)
  }

  initial <- .gtMaxwellModes(
    log(tau), t, Gexp, wexp, isPlateau
  )
  g <- initial$g
  success <- FALSE

  # A plateau-only solution has no relaxation times to optimize.
  if (isTRUE(isPlateau) && length(tau) == 0L) {
    return(list(success = TRUE, g = g, tau = tau))
  }

  initErr <- sqrt(sum(res_tG(log(tau))^2))

  if (requireNamespace("minpack.lm", quietly = TRUE)) {
    tryCatch({
      fit <- minpack.lm::nls.lm(
        par = log(tau),
        fn = res_tG,
        control = minpack.lm::nls.lm.control(maxiter = 200)
      )

      tau_new <- exp(fit$par)

      res_new <- .gtMaxwellModes(
        log(tau_new), t, Gexp, wexp, isPlateau
      )

      finalErr <- sqrt(sum(
        res_tG(log(res_new$tau))^2
      ))

      if (finalErr <= initErr) {
        tau <- res_new$tau
        g <- res_new$g
        success <- TRUE
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
    isPlateau <- isTRUE(par$plateau)

  if (isPlateau) {
    if (is.null(crs$G0) ||
        length(crs$G0) != 1L ||
        !is.finite(crs$G0)) {
      stop(
        "Time-domain plateau fitting requires a continuous result with G0.",
        call. = FALSE
      )
    }
  }

  Nmin <- max(floor(0.5 * log10(max(t) / min(t))), 2)
  Nmax <- min(floor(3.0 * log10(max(t) / min(t))), n / 4)
  if (par$MaxNumModes > 0) Nmax <- min(Nmax, par$MaxNumModes)
  Nv   <- as.integer(seq(Nmin, Nmax))

    Gc <- .gtKernelPrestore(
    H, kernMat,
    if (isPlateau) crs$G0 else NULL
  )
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
      res <- .gtMaxwellModes(
        gd$z, t, Gexp, wexp, isPlateau
      )
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
  res <- .gtMaxwellModes(
        gd$z, t, Gexp, wexp, isPlateau
      )
  g   <- res$g; tau <- res$tau; error <- res$error; cKp <- res$condKp

    ft <- .gtFineTuneSolution(
    tau, t, Gexp, wexp, isPlateau
  )
  if (ft$success) { g <- ft$g; tau <- ft$tau }

    idx <- order(tau)
  tau <- tau[idx]

  if (isPlateau) {
    plateau_weight <- g[length(g)]
    g <- c(g[idx], plateau_weight)
  } else {
    g <- g[idx]
  }
  if (length(tau) > 1) {
    tauSpacing <- tau[-1] / tau[-length(tau)]
    itry <- 0
    while (min(tauSpacing) < par$minTauSpacing && itry < 3) {
      imode <- which.min(tauSpacing)
            if (isPlateau) {
        plateau_weight <- g[length(g)]

        mg <- .gtMergeModes(
          g[seq_along(tau)], tau, imode
        )

        g <- c(mg$g, plateau_weight)
        tau <- mg$tau
      } else {
        mg <- .gtMergeModes(g, tau, imode)
        g <- mg$g
        tau <- mg$tau
      }
        ft <- .gtFineTuneSolution(
		tau, t, Gexp, wexp, isPlateau
		)
      if (ft$success) { g <- ft$g; tau <- ft$tau }
      if (length(tau) > 1) tauSpacing <- tau[-1] / tau[-length(tau)] else break
      itry <- itry + 1
    }
  }

  # Report the weighted SSE of the final returned modes.
  G0 <- NULL

  if (isPlateau) {
    stopifnot(length(g) == length(tau) + 1L)
    G0 <- g[length(g)]
    g <- g[seq_along(tau)]
  }

  # Report the weighted SSE of the final returned modes and plateau.
  Kfinal <- exp(-outer(t, 1 / tau))
  Gfinal <- as.vector(Kfinal %*% g)

  if (isPlateau) Gfinal <- Gfinal + G0

  error <- sum((wexp * (Gfinal / Gexp - 1))^2)

  if (isTRUE(writeOutput)) {
    if (!dir.exists(outputDir)) {
      dir.create(outputDir, recursive = TRUE)
    }

    modes_file <- file.path(outputDir, "dmodes.dat")

    if (isPlateau) {
      write(sprintf("# G0 = %.17e", G0), modes_file)
      utils::write.table(
        cbind(g, tau), modes_file,
        append = TRUE,
        row.names = FALSE,
        col.names = FALSE
      )
    } else {
      utils::write.table(
        cbind(g, tau), modes_file,
        row.names = FALSE,
        col.names = FALSE
      )
    }

    utils::write.table(
      cbind(wtBase, nzNbst, AICbst),
      file.path(outputDir, "aic.dat"),
      row.names = FALSE,
      col.names = FALSE
    )

    utils::write.table(
      cbind(t, Gfinal),
      file.path(outputDir, "Gfitd.dat"),
      row.names = FALSE,
      col.names = FALSE
    )
  }

    list(
    g = g,
    tau = tau,
    G0 = G0,
    error = error,
    condKp = cKp,
    Nopt = Nopt
  )

}
