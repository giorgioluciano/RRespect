
# MaxwellModes.R


library(nnls)
library(pracma)

nnLLS <- function(w, tau, Gexp, wexp) {
  n    <- length(w)
  ntau <- length(tau)

  res  <- meshgrid(tau, w)   # X = tau (cols), Y = w (rows)  → ws = X*Y
  S    <- res$X; W <- res$Y
  ws   <- S * W; ws2 <- ws^2

  K  <- rbind(ws2 / (1 + ws2), ws / (1 + ws2))   # 2n x ntau

  # weighted: Kp_(i,j) = (wexp_i / Gexp_i) * K_(i,j)
  Kp <- (wexp / Gexp) * K
  condKp <- pracma::cond(Kp)

  g_sol  <- nnls::nnls(Kp, wexp)$x

  GstM  <- K %*% g_sol
  error <- sum((wexp * (as.vector(GstM) / Gexp - 1))^2)

  return(list(g = g_sol, error = error, condKp = condKp))
}

MaxwellModes <- function(z, w, Gexp, wexp) {
  tau <- exp(z)
  n   <- length(w)

  res    <- nnLLS(w, tau, Gexp, wexp)
  g      <- res$g
  error  <- res$error
  condKp <- res$condKp

  # remove runaway modes outside frequency window
  izero <- which(max(w) * min(tau) < 0.02 | min(w) * max(tau) > 50)
  if (length(izero) > 0) { tau <- tau[-izero]; g <- g[-izero] }

  # prune negligible weights
  if (length(g) > 0 && max(g) > 0) {
    izero <- which(g / max(g) < 1e-8)
    if (length(izero) > 0) { tau <- tau[-izero]; g <- g[-izero] }
  }

  return(list(g = g, tau = tau, error = error, condKp = condKp))
}
