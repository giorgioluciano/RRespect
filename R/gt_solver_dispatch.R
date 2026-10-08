
# Continuous time-domain solver dispatch.
.gtFitH <- function(
  lam, Gexp, wexp, H, kernMat, G0 = NULL,
  solver = "legacy", max_nfev = NULL
) {
  solver <- match.arg(solver, c("legacy", "experimental"))

  if (identical(solver, "legacy")) {
    out <- .gtLevenMarq(lam, Gexp, wexp, H, kernMat, G0)
    return(list(
      H = if (is.null(G0)) out else out$H,
      G0 = if (is.null(G0)) NULL else out$G0,
      diagnostics = NULL,
      solver = "legacy"
    ))
  }

  out <- .gtGetHExperimental(
    lam = lam,
    Gexp = Gexp,
    wexp = wexp,
    H = H,
    kernMat = kernMat,
    G0 = G0,
    max_nfev = max_nfev
  )

  if (!isTRUE(out$diagnostics$success)) {
    stop(
      sprintf(
        paste0(
          "Temporal TRF did not converge at lambda = %.6e ",
          "(status = %d, nfev = %d)."
        ),
        lam,
        out$diagnostics$status,
        out$diagnostics$nfev
      ),
      call. = FALSE
    )
  }

  out
}

.gtInitializeFit <- function(
  Gexp, wexp, s, kernMat, G0 = NULL, solver = "legacy"
) {
  H <- -5 * rep(1, length(s)) + sin(pi * s)
  .gtFitH(
    lam = 1,
    Gexp = Gexp,
    wexp = wexp,
    H = H,
    kernMat = kernMat,
    G0 = G0,
    solver = solver
  )
}
