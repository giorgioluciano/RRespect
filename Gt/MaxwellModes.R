
# MaxwellModes_time.R
# Time-domain DRS: kernel K_(i,j) = exp(-t_i / tau_j)

library(nnls)
library(pracma)
source("common.R")

nnLLS_time <- function(t, tau, Gexp, wexp) {
  res  <- meshgrid(tau, t)
  S    <- res$X;  T <- res$Y
  K    <- exp(-T / S)                              # n x ntau

  Kp     <- (wexp / Gexp) * K                     # weighted rows
  condKp <- pracma::cond(Kp)
  g_sol  <- nnls::nnls(Kp, wexp)$x

  GtM   <- K %*% g_sol
  error <- sum((wexp * (as.vector(GtM) / Gexp - 1))^2)
  return(list(g = g_sol, error = error, condKp = condKp))
}

MaxwellModes <- function(z, t, Gexp, wexp) {
  tau <- exp(z)
  res    <- nnLLS_time(t, tau, Gexp, wexp)
  g      <- res$g;  error <- res$error;  condKp <- res$condKp

  # prune negligible weights
  if (length(g) > 0 && max(g) > 0) {
    izero <- which(g / max(g) < 1e-7)
    if (length(izero) > 0) { tau <- tau[-izero]; g <- g[-izero] }
  }
  return(list(g = g, tau = tau, error = error, condKp = condKp))
}
