# gt_levenmarq.R — Levenberg-Marquardt solver, time domain (internal)
# Ported from Gt/LevenMarq.R; no library() / source() calls.

.gtKernelD <- function(H, kernMat) {
  n  <- nrow(kernMat)
  ns <- ncol(kernMat)
  Hsuper <- matrix(exp(H), nrow = n, ncol = ns, byrow = TRUE)
  kernMat * Hsuper
}

.gtResidualLM <- function(H, lam, Gexp, wexp, kernMat) {
  n  <- nrow(kernMat)
  ns <- ncol(kernMat)
  nl <- ns - 2

  G0 <- NULL
  if (length(H) == ns + 1L) {
    G0 <- H[ns + 1L]
    H <- H[seq_len(ns)]
  }

  r <- numeric(n + nl)
  r[seq_len(n)] <- wexp * (
    1 - .gtKernelPrestore(H, kernMat, G0) / Gexp
  )
  r[n + seq_len(nl)] <- sqrt(lam) * diff(H, differences = 2)
  r
}
.gtJacobianLM <- function(H, lam, Gexp, wexp, kernMat) {
  n  <- nrow(kernMat)
  ns <- ncol(kernMat)
  nl <- ns - 2

  plateau <- length(H) == ns + 1L
  if (plateau) H <- H[seq_len(ns)]

  idx <- seq_len(nl)
  L <- matrix(0, nl, ns)
  L[cbind(idx, idx)]     <-  1
  L[cbind(idx, idx + 1)] <- -2
  L[cbind(idx, idx + 2)] <-  1

  Kmatrix <- matrix(wexp / Gexp, nrow = n, ncol = ns)
  Jr <- matrix(0, n + nl, ns + as.integer(plateau))

  Jr[seq_len(n), seq_len(ns)] <-
    -.gtKernelD(H, kernMat) * Kmatrix

  Jr[n + seq_len(nl), seq_len(ns)] <- sqrt(lam) * L

  if (plateau) {
    Jr[seq_len(n), ns + 1L] <- -wexp / Gexp
  }

  Jr
}

.gtLevenMarq <- function(lam, Gexp, wexp, H, kernMat, G0 = NULL) {
  ns <- ncol(kernMat)
  plateau <- !is.null(G0)

  if (plateau) H <- c(H, G0)
  npar <- length(H)

  nu   <- 2
  iter <- 0

  r   <- .gtResidualLM(H, lam, Gexp, wexp, kernMat)
  Jr  <- .gtJacobianLM(H, lam, Gexp, wexp, kernMat)
  JtJ <- crossprod(Jr)
  mu  <- 1e-3 * max(diag(JtJ))

  repeat {
    iter <- iter + 1

    A <- JtJ + mu * diag(npar)
    rhs   <- -as.vector(t(Jr) %*% r)
    Delta <- tryCatch(solve(A, rhs),error = function(e) rep(0, npar))

    Hnew  <- H + Delta
    r_new <- .gtResidualLM(Hnew, lam, Gexp, wexp, kernMat)

    rho_num <- sum(r^2) - sum(r_new^2)
    rho_den <- as.numeric(t(Delta) %*% (mu * Delta + rhs))

    rho_ok  <- is.finite(rho_den) && abs(rho_den) > 1e-15
    rho_val <- if (rho_ok) rho_num / rho_den else -1

    if (is.finite(rho_val) && rho_val > 0) {
      H   <- Hnew
      r   <- r_new
      Jr  <- .gtJacobianLM(H, lam, Gexp, wexp, kernMat)
      JtJ <- crossprod(Jr)
      mu  <- mu * max(1/3, 1 - (2 * rho_val - 1)^3)
      nu  <- 2

      grad_norm  <- max(abs(t(Jr) %*% r))
      delta_norm <- sqrt(sum(Delta^2))
      H_norm     <- sqrt(sum(H^2))

      if (grad_norm < 1e-6 || delta_norm < 1e-6 * (1 + H_norm)) break
    } else {
      mu <- mu * nu
      nu <- 2 * nu
      if (mu > 1e15) break
    }

    if (iter >= 5000) break
  }
    if (plateau) {
    return(list(
      H = H[seq_len(ns)],
      G0 = H[ns + 1L]
    ))
  }

  H
}
