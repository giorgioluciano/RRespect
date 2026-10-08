
.gtGetHExperimental <- function(
  lam, Gexp, wexp, H, kernMat, G0 = NULL,
  ftol = 1e-8, xtol = 1e-8, gtol = 1e-8,
  max_nfev = NULL
) {
  stopifnot(
    length(lam) == 1L,
    is.finite(lam),
    lam >= 0,
    is.matrix(kernMat),
    is.numeric(kernMat),
    nrow(kernMat) >= 1L,
    ncol(kernMat) >= 3L,
    all(is.finite(kernMat)),
    length(Gexp) == nrow(kernMat),
    all(is.finite(Gexp)),
    all(Gexp > 0),
    length(wexp) == length(Gexp),
    all(is.finite(wexp)),
    all(wexp >= 0),
    any(wexp > 0),
    length(H) == ncol(kernMat),
    all(is.finite(H))
  )

  if (!is.null(G0)) {
    stopifnot(length(G0) == 1L, is.finite(G0))
  }

  ns <- ncol(kernMat)
  start <- if (is.null(G0)) H else c(H, G0)

  fit <- .trfNoBounds(
    fun = function(x) {
      .gtResidualLM(x, lam, Gexp, wexp, kernMat)
    },
    jac = function(x) {
      .gtJacobianLM(x, lam, Gexp, wexp, kernMat)
    },
    x0 = start,
    ftol = ftol,
    xtol = xtol,
    gtol = gtol,
    max_nfev = max_nfev
  )

  Hopt <- fit$x[seq_len(ns)]
  G0opt <- if (is.null(G0)) NULL else fit$x[ns + 1L]

  list(
    H = Hopt,
    G0 = G0opt,
    lam = lam,
    Gfit = .gtKernelPrestore(Hopt, kernMat, G0opt),
    diagnostics = fit,
    solver = "experimental"
  )
}
