
# LevenMarq_time.R
library(pracma)
source("common.R")

kernelD_time <- function(H, kernMat) {
  n  <- nrow(kernMat)
  ns <- ncol(kernMat)
  Hsuper <- matrix(exp(H), nrow = n, ncol = ns, byrow = TRUE)
  return(kernMat * Hsuper)
}

residualLM_time <- function(H, lam, Gexp, wexp, kernMat) {
  n  <- nrow(kernMat)
  ns <- ncol(kernMat)
  nl <- ns - 2
  r  <- numeric(n + nl)
  r[1:n]          <- wexp * (1 - kernel_prestore(H, kernMat) / Gexp)
  r[(n+1):(n+nl)] <- sqrt(lam) * diff(H, differences = 2)
  return(r)
}

jacobianLM_time <- function(H, lam, Gexp, wexp, kernMat) {
  n  <- nrow(kernMat)
  ns <- ncol(kernMat)
  nl <- ns - 2

  L <- matrix(0, nl, ns)
  for (i in 1:nl) { L[i,i] <- 1; L[i,i+1] <- -2; L[i,i+2] <- 1 }

  Kmatrix <- matrix(wexp / Gexp, nrow = n, ncol = ns)
  Jr <- matrix(0, n + nl, ns)
  Jr[1:n,          1:ns] <- -kernelD_time(H, kernMat) * Kmatrix
  Jr[(n+1):(n+nl), 1:ns] <- sqrt(lam) * L
  return(Jr)
}

LevenMarq <- function(lam, Gexp, wexp, H, kernMat) {
  ns   <- ncol(kernMat)
  nu   <- 2
  iter <- 0

  r   <- residualLM_time(H, lam, Gexp, wexp, kernMat)
  Jr  <- jacobianLM_time(H, lam, Gexp, wexp, kernMat)
  JtJ <- crossprod(Jr)
  mu  <- 1e-3 * max(diag(JtJ))

  repeat {
    iter <- iter + 1

    # solve (JtJ + mu*I) * Delta = -Jt*r
    A     <- JtJ + mu * diag(ns)
    rhs   <- -as.vector(t(Jr) %*% r)
    Delta <- tryCatch(solve(A, rhs), error = function(e) rep(0, ns))

    Hnew  <- H + Delta
    r_new <- residualLM_time(Hnew, lam, Gexp, wexp, kernMat)

    rho_num <- sum(r^2) - sum(r_new^2)
    rho_den <- as.numeric(t(Delta) %*% (mu * Delta + rhs))

    # robusto: se den e' NA/zero/NaN tratta come rho < 0
    rho_ok  <- is.finite(rho_den) && abs(rho_den) > 1e-15
    rho_val <- if (rho_ok) rho_num / rho_den else -1

    if (is.finite(rho_val) && rho_val > 0) {
      H   <- Hnew
      r   <- r_new
      Jr  <- jacobianLM_time(H, lam, Gexp, wexp, kernMat)
      JtJ <- crossprod(Jr)
      mu  <- mu * max(1/3, 1 - (2*rho_val - 1)^3)
      nu  <- 2

      grad_norm  <- max(abs(t(Jr) %*% r))
      delta_norm <- sqrt(sum(Delta^2))
      H_norm     <- sqrt(sum(H^2))

      if (grad_norm < 1e-6 || delta_norm < 1e-6 * (1 + H_norm)) break
    } else {
      mu <- mu * nu
      nu <- 2 * nu
      if (mu > 1e15) break   # evita loop infinito se diverge
    }

    if (iter >= 5000) break
  }
  return(H)
}
