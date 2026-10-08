
# Python v2.1 error-band calculation for experimental workflows.
.respectErrorBand <- function(H_lam, p_lam) {
  stopifnot(
    is.matrix(H_lam),
    is.numeric(H_lam),
    nrow(H_lam) >= 1L,
    ncol(H_lam) >= 1L,
    all(is.finite(H_lam)),
    is.numeric(p_lam),
    length(p_lam) == ncol(H_lam),
    all(is.finite(p_lam)),
    all(p_lam >= 0),
    abs(sum(p_lam) - 1) < 1e-8
  )

  selected <- which(p_lam > 0.1)

  if (length(selected)) {
    Hm <- numeric(nrow(H_lam))
    Hm2 <- numeric(nrow(H_lam))

    for (i in selected) {
      Hm <- Hm + H_lam[, i]
      Hm2 <- Hm2 + H_lam[, i]^2
    }

    Hm <- Hm / length(selected)
    Hm2 <- Hm2 / length(selected)
  } else {
    Hm <- as.vector(H_lam %*% p_lam)
    Hm2 <- as.vector((H_lam^2) %*% p_lam)
  }

  sqrt(pmax(Hm2 - Hm^2, 0))
}
