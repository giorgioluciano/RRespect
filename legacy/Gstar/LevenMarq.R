# LevenMarq.R
# Ottimizzazione Levenberg-Marquardt per H (e opzionalmente G0)

levenberg_marquardt <- function(par, s, Gstar, kernMat, H0, lam, G0 = NULL) {

  n      <- length(Gstar) / 2        # numero di frequenze
  ns     <- length(H0)
  H      <- H0
  use_G0 <- !is.null(G0)

  ds <- diff(log(s))                 # pesi integrazione trapezi (log-s)
  ds <- c(ds[1], ds, ds[length(ds)])
  ds <- 0.5 * (ds[-length(ds)] + ds[-1])

  # Funzione residuo r(H) = K*H - Gstar
  residual <- function(H, G0 = NULL) {
    r <- kernMat %*% H - Gstar
    if (!is.null(G0)) r[1:n] <- r[1:n] + G0
    r
  }

  # Jacobiano: K (+ colonna 1 per G0)
  J_full <- if (use_G0) cbind(kernMat, c(rep(1, n), rep(0, n))) else kernMat

  H_aug  <- if (use_G0) c(H, G0) else H
  ns_aug <- length(H_aug)

  mu     <- lam
  nu     <- 2
  iter   <- 0
  converged <- FALSE

  while (iter < par$MaxIter && !converged) {
    iter  <- iter + 1
    r     <- if (use_G0) residual(H_aug[1:ns], H_aug[ns + 1]) else residual(H_aug)
    JtJ   <- crossprod(J_full)          # J'J
    Jtr   <- crossprod(J_full, r)       # J'r

    # Regolarizzazione L2 su H (non su G0)
    reg            <- diag(ns_aug) * mu
    if (use_G0) reg[ns_aug, ns_aug] <- 0   # non regolarizza G0

    delta <- tryCatch(
      solve(JtJ + reg, -Jtr),
      error = function(e) rep(0, ns_aug)
    )

    H_new <- H_aug + delta
    H_new[1:ns] <- pmax(H_new[1:ns], 0)   # H >= 0

    r_new <- if (use_G0) residual(H_new[1:ns], H_new[ns + 1]) else residual(H_new)
    rho   <- (sum(r^2) - sum(r_new^2)) /
             (sum(delta * (mu * delta - Jtr)))

    if (rho > 0) {
      H_aug <- H_new
      mu    <- mu * max(1/3, 1 - (2*rho - 1)^3)
      nu    <- 2
    } else {
      mu <- mu * nu
      nu <- nu * 2
    }

    if (max(abs(delta)) < par$Ptol || sum(r^2) < par$Ftol) converged <- TRUE
  }

  out <- list(H = H_aug[1:ns])
  if (use_G0) out$G0 <- H_aug[ns + 1]
  return(out)
}
