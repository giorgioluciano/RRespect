# gt_maxwellmodes.R — Maxwell mode fitting, time domain (internal)
# Ported from Gt/MaxwellModes.R; uses .gridDensity from gridDensity_shared.R.

.gtNnLLS <- function(t, tau, Gexp, wexp, isPlateau = FALSE) {
  K <- exp(-outer(t, 1 / tau))

  if (isTRUE(isPlateau)) {
    K <- cbind(K, rep(1, length(t)))
  }

  Kp     <- (wexp / Gexp) * K
  condKp <- pracma::cond(Kp)
  g_sol  <- nnls::nnls(Kp, wexp)$x

  GtM   <- as.vector(K %*% g_sol)
  error <- sum((wexp * (GtM / Gexp - 1))^2)

  list(
    g = g_sol,
    error = error,
    condKp = condKp
  )
}

.gtMaxwellModes <- function(z, t, Gexp, wexp, isPlateau = FALSE) {
  tau <- exp(z)

  res <- .gtNnLLS(t, tau, Gexp, wexp, isPlateau)
  g <- res$g
  error <- res$error
  condKp <- res$condKp

  if (isTRUE(isPlateau)) {
    g_main <- g[seq_along(tau)]
  } else {
    g_main <- g
  }

  if (length(g_main) > 0 && max(g_main) > 0) {
    izero <- which(g_main / max(g_main) < 1e-7)

    if (length(izero) > 0) {
      tau <- tau[-izero]
      g <- g[-izero]
    }
  }

  list(
    g = g,
    tau = tau,
    error = error,
    condKp = condKp
  )
}