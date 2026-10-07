.gridDensity <- function(x, px, N) {
  npts <- 100
  xi <- seq(min(x), max(x), length.out = npts)
  pint <- pracma::interp1(as.numeric(x), as.numeric(px), xi, method = "spline")
  pint[is.na(pint) | pint < 0] <- 0

  ci <- pracma::cumtrapz(xi, pint)
  total <- ci[npts]
  if (total <= 0) total <- 1
  ci <- ci / total
  ci[1] <- 0
  ci[npts] <- 1

  keep <- c(TRUE, diff(ci) > 1e-12)
  ci_k <- ci[keep]
  xi_k <- xi[keep]

  fint_inv <- function(b) {
    b <- pmax(pmin(b, max(ci_k)), min(ci_k))
    pracma::interp1(as.numeric(ci_k), as.numeric(xi_k), as.numeric(b), method = "spline")
  }

  alfa <- 1 / (N - 1)
  z <- numeric(N)
  zij <- numeric(N + 1)

  z[1] <- min(x)
  z[N] <- max(x)
  zij[1] <- z[1]
  zij[N + 1] <- z[N]

  if (N > 2) {
    beta_ij <- seq(0.5, N - 1.5) * alfa
    zij[2:N] <- vapply(beta_ij, fint_inv, numeric(1))
    beta_z <- seq(1, N - 2) * alfa
    z[2:(N - 1)] <- vapply(beta_z, fint_inv, numeric(1))
  }

  h <- diff(zij)
  list(z = z, h = h)
}
