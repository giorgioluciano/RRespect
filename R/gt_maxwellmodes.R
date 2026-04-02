# gt_maxwellmodes.R — Maxwell mode fitting, time domain (internal)
# Ported from Gt/MaxwellModes.R; uses .gridDensity from gridDensity_shared.R.

.gtNnLLS <- function(t, tau, Gexp, wexp) {
  K    <- exp(-outer(t, 1/tau))                # n x ntau: K[i,j] = exp(-t[i]/tau[j])

  Kp     <- (wexp / Gexp) * K                 # weighted rows
  condKp <- pracma::cond(Kp)
  g_sol  <- nnls::nnls(Kp, wexp)$x

  GtM   <- K %*% g_sol
  error <- sum((wexp * (as.vector(GtM) / Gexp - 1))^2)
  list(g = g_sol, error = error, condKp = condKp)
}

.gtMaxwellModes <- function(z, t, Gexp, wexp) {
  tau <- exp(z)
  res <- .gtNnLLS(t, tau, Gexp, wexp)
  g   <- res$g;  error <- res$error;  condKp <- res$condKp

  if (length(g) > 0 && max(g) > 0) {
    izero <- which(g / max(g) < 1e-7)
    if (length(izero) > 0) { tau <- tau[-izero]; g <- g[-izero] }
  }
  list(g = g, tau = tau, error = error, condKp = condKp)
}
