
.experimentalPseudoInverse <- function(A, rcond = 1e-15) {
  decomposition <- svd(A)
  keep <- decomposition$d > rcond * max(decomposition$d)

  if (!any(keep)) {
    return(matrix(0, ncol(A), nrow(A)))
  }

  sweep(
    decomposition$v[, keep, drop = FALSE],
    2, decomposition$d[keep], "/"
  ) %*% t(decomposition$u[, keep, drop = FALSE])
}

.experimentalFineTune <- function(
  tau, axis, Gexp, wexp, plateau, domain
) {
  residual <- function(candidate) {
    fit <- .experimentalNnLS(
      candidate, axis, Gexp, wexp, plateau, domain
    )
    wexp * (fit$fitted / Gexp - 1)
  }

  initial_error <- sqrt(sum(residual(tau)^2))
  bounds <- .experimentalTauBounds(axis, domain)

  attempt <- tryCatch({
    fit <- .trfBounds(
      fun = residual,
      jac = function(candidate) {
        .experimentalJacobian2Point(
          residual, candidate, bounds[1L], bounds[2L]
        )
      },
      x0 = tau,
      lb = bounds[1L],
      ub = bounds[2L]
    )

    covariance <- .experimentalPseudoInverse(
      crossprod(fit$jac)
    ) * mean(fit$fun^2)

    dtau <- sqrt(diag(covariance))

    list(fit = fit, dtau = dtau)
  }, error = function(e) e)

  if (inherits(attempt, "error")) {
    warning(
      paste0(
        "Discrete TRF fine-tuning failed; keeping the NNLS solution: ",
        conditionMessage(attempt)
      ),
      call. = FALSE
    )

    starting <- .experimentalMaxwellModes(
      log(tau), axis, Gexp, wexp, plateau, domain
    )

    return(list(
      accepted = FALSE,
      g = starting$g,
      tau = starting$tau,
      dtau = rep(NA_real_, length(starting$tau)),
      diagnostics = NULL
    ))
  }

  modes <- .experimentalMaxwellModes(
    log(attempt$fit$x), axis, Gexp, wexp, plateau, domain
  )

  final_error <- sqrt(sum(residual(modes$tau)^2))

  list(
    accepted = final_error <= initial_error,
    g = modes$g,
    tau = modes$tau,
    dtau = attempt$dtau[modes$keep],
    diagnostics = attempt$fit
  )
}

.experimentalMergeModes <- function(g, tau, imode, domain) {
  g1 <- g[imode]
  g2 <- g[imode + 1L]
  tau1 <- tau[imode]
  tau2 <- tau[imode + 1L]

  cost <- function(parameters) {
    gn <- parameters[1L]
    taun <- parameters[2L]

    if (identical(domain, "time")) {
      integrand_time <- function(t) {
        new <- gn * exp(-t / taun)
        old <- g1 * exp(-t / tau1) + g2 * exp(-t / tau2)
        (new / old - 1)^2
      }

      lower <- min(tau1, tau2) / 10
      upper <- max(tau1, tau2) * 10
    } else {
      integrand_frequency <- function(w) {
        wt <- w * taun
        new_p <- gn * wt^2 / (1 + wt^2)
        new_pp <- gn * wt / (1 + wt^2)

        wt1 <- w * tau1
        wt2 <- w * tau2

        old_p <- g1 * wt1^2 / (1 + wt1^2) +
          g2 * wt2^2 / (1 + wt2^2)
        old_pp <- g1 * wt1 / (1 + wt1^2) +
          g2 * wt2 / (1 + wt2^2)

        (new_p / old_p - 1)^2 + (new_pp / old_pp - 1)^2
      }

      lower <- min(1 / tau1, 1 / tau2) / 10
      upper <- max(1 / tau1, 1 / tau2) * 10
    }

    stats::integrate(if (identical(domain, "time")) integrand_time else integrand_frequency, lower, upper,
      subdivisions = 50L,
      rel.tol = 1.49e-8,
      abs.tol = 1.49e-8
    )$value
  }

  gradient <- function(parameters) {
    f0 <- cost(parameters)
    h <- sqrt(.Machine$double.eps)

    vapply(seq_along(parameters), function(j) {
      perturbed <- parameters
      perturbed[j] <- parameters[j] + h
      (cost(perturbed) - f0) /
        (perturbed[j] - parameters[j])
    }, numeric(1))
  }

  fit <- stats::optim(
    c(g1 + g2, 0.5 * (tau1 + tau2)),
    fn = cost,
    gr = gradient,
    method = "BFGS",
    control = list(maxit = 400L, reltol = 1e-12)
  )

  merged <- tau[-(imode + 1L)]
  merged[imode] <- fit$par[2L]
  merged
}

.experimentalModeCounts <- function(axis, max_modes) {
  decades <- log10(max(axis) / min(axis))
  minimum <- as.integer(max(floor(0.5 * decades), 2))
  maximum <- as.integer(min(floor(3 * decades), length(axis) / 4))

  if (!is.null(max_modes) && max_modes > 0) {
    maximum <- min(maximum, as.integer(max_modes))
  }

  if (maximum < 1L) {
    stop("The data window does not support a discrete mode scan.")
  }

  if (maximum > minimum) {
    return(seq.int(minimum, maximum))
  }

  warning(
    sprintf(
      "Only N = %d is scanned; check maxNumModes and the data window.",
      maximum
    ),
    call. = FALSE
  )

  maximum
}

.experimentalDiscSpec <- function(
  par, crs, domain, writeOutput = FALSE, outputDir = "output"
) {
  if (!identical(crs$solver, "experimental")) {
    stop(
      "Experimental discrete fitting requires an experimental continuous result.",
      call. = FALSE
    )
  }

  axis <- if (identical(domain, "time")) crs$t else crs$w
  Gexp <- crs$Gexp
  wexp <- crs$wexp
  plateau <- isTRUE(par$plateau)
  s <- crs$s
  H <- crs$H

  Nv <- .experimentalModeCounts(axis, par$MaxNumModes)
  Cerror <- 1 / .experimentalPopulationSD(
    wexp * (crs$G_fit / Gexp - 1)
  )

  if (!is.finite(Cerror)) {
    stop("Non-finite discrete AIC error weight.", call. = FALSE)
  }

  delta <- par$deltaBaseWeightDist
  base_count <- as.integer(1 / delta) - 1L

  if (base_count < 1L) {
    stop("No base weights are available for the discrete AIC scan.")
  }

  wt_base <- delta * seq_len(base_count)
  AIC_bst <- numeric(length(wt_base))
  N_bst <- nz_N_bst <- integer(length(wt_base))

  for (ib in seq_along(wt_base)) {
    wt <- .experimentalPlacementWeights(
      H, axis, s, wt_base[ib], domain
    )

    errors <- numeric(length(Nv))
    surviving <- integer(length(Nv))

    for (i in seq_along(Nv)) {
      grid <- .experimentalGridDensity(log(s), wt, Nv[i])
      modes <- .experimentalMaxwellModes(
        grid$z, axis, Gexp, wexp, plateau, domain
      )

      errors[i] <- modes$error
      surviving[i] <- length(modes$tau)
    }

    AIC <- 2 * Nv + 2 * Cerror * errors
    best <- which.min(AIC)

    AIC_bst[ib] <- AIC[best]
    N_bst[ib] <- Nv[best]
    nz_N_bst[ib] <- surviving[best]
  }

  best <- which.min(AIC_bst)
  wt <- .experimentalPlacementWeights(
    H, axis, s, wt_base[best], domain
  )
  grid <- .experimentalGridDensity(log(s), wt, N_bst[best])
  modes <- .experimentalMaxwellModes(
    grid$z, axis, Gexp, wexp, plateau, domain
  )

  g <- modes$g
  tau <- modes$tau
  dtau <- rep(NA_real_, length(tau))
  refinements <- list()

  refinement <- .experimentalFineTune(
    tau, axis, Gexp, wexp, plateau, domain
  )
  refinements[[1L]] <- refinement$diagnostics

  if (isTRUE(refinement$accepted)) {
    g <- refinement$g
    tau <- refinement$tau
    dtau <- refinement$dtau
  }

  merges <- 0L

  while (length(tau) > 1L && merges < 3L) {
    spacing <- tau[-1L] / tau[-length(tau)]
    if (min(spacing) >= par$minTauSpacing) break

    merged_tau <- .experimentalMergeModes(
      g, tau, which.min(spacing), domain
    )

    refinement <- .experimentalFineTune(
      merged_tau, axis, Gexp, wexp, plateau, domain
    )
    refinements[length(refinements) + 1L] <-
      list(refinement$diagnostics)

    if (isTRUE(refinement$accepted)) {
      g <- refinement$g
      tau <- refinement$tau
      dtau <- refinement$dtau
    } else {
      modes <- .experimentalMaxwellModes(
        log(merged_tau), axis, Gexp, wexp, plateau, domain
      )
      g <- modes$g
      tau <- modes$tau
      dtau <- rep(NA_real_, length(tau))
    }

    merges <- merges + 1L
  }

  G0 <- 0
  if (plateau) {
    G0 <- utils::tail(g, 1L)
    g <- utils::head(g, -1L)
  }

  G_fit <- as.vector(
    .experimentalMaxwellKernel(tau, axis, domain) %*% g
  )

  if (identical(domain, "time")) {
    G_fit <- G_fit + G0
  } else {
    G_fit[seq_along(axis)] <- G_fit[seq_along(axis)] + G0
  }

  error <- sum((wexp * (G_fit / Gexp - 1))^2)

  out <- list(
    g = g,
    tau = tau,
    dtau = dtau,
    N = length(g),
    Nopt = N_bst[best],
    G0 = G0,
    G_fit = G_fit,
    error = error,
    wt_base = wt_base,
    AIC_bst = AIC_bst,
    N_bst = N_bst,
    nz_N_bst = nz_N_bst,
    solver = "experimental",
    reference = crs$reference,
    diagnostics = list(
      Cerror = Cerror,
      mode_counts = Nv,
      merges = merges,
      refinements = refinements
    )
  )

  out$Gfit <- if (identical(domain, "time")) {
    cbind(axis, G_fit)
  } else {
    n <- length(axis)
    cbind(axis, G_fit[seq_len(n)], G_fit[n + seq_len(n)])
  }

  if (isTRUE(writeOutput)) {
    dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)

    header <- if (plateau) {
      c(sprintf("G0 = %.6e", G0), "g  tau  dtau")
    } else {
      "g  tau  dtau"
    }

    file <- file.path(outputDir, "drs.dat")
    writeLines(paste0("# ", header), file)
    utils::write.table(
      cbind(g, tau, dtau),
      file,
      append = TRUE,
      row.names = FALSE,
      col.names = FALSE
    )

    file <- file.path(outputDir, "aic.dat")
    writeLines("# wt_base  N_bst  AIC  nz_N_bst", file)
    utils::write.table(
      cbind(wt_base, N_bst, AIC_bst, nz_N_bst),
      file,
      append = TRUE,
      row.names = FALSE,
      col.names = FALSE
    )
  }

  out
}
