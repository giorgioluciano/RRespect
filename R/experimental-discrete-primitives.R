
.experimentalMaxwellKernel <- function(tau, axis, domain) {
  if (identical(domain, "time")) {
    return(exp(-outer(axis, 1 / tau)))
  }

  ws <- outer(axis, tau)
  rbind(ws^2 / (1 + ws^2), ws / (1 + ws^2))
}

.experimentalTauBounds <- function(axis, domain) {
  if (identical(domain, "time")) {
    c(0.02 * min(axis), 50 * max(axis))
  } else {
    c(0.02 / max(axis), 50 / min(axis))
  }
}

.experimentalNnLS <- function(
  tau, axis, Gexp, wexp, plateau, domain
) {
  K <- .experimentalMaxwellKernel(tau, axis, domain)

  if (plateau) {
    column <- if (identical(domain, "time")) {
      rep(1, length(axis))
    } else {
      c(rep(1, length(axis)), rep(0, length(axis)))
    }
    K <- cbind(K, column)
  }

  stopifnot(ncol(K) > 0L)
  weighted_K <- sweep(K, 1, wexp / Gexp, "*")
  g <- as.vector(nnls::nnls(weighted_K, wexp)$x)
  fitted <- as.vector(K %*% g)

  list(
    g = g,
    error = sum((wexp * (fitted / Gexp - 1))^2),
    fitted = fitted
  )
}

.experimentalMaxwellModes <- function(
  z, axis, Gexp, wexp, plateau, domain
) {
  tau <- exp(z)
  fit <- .experimentalNnLS(tau, axis, Gexp, wexp, plateau, domain)
  g_modes <- fit$g[seq_along(tau)]
  bounds <- .experimentalTauBounds(axis, domain)

  keep <- which(tau >= bounds[1L] & tau <= bounds[2L])

  if (length(keep)) {
    reference <- max(g_modes[keep])
    if (reference > 0) {
      keep <- keep[g_modes[keep] / reference >= 1e-7]
    } else {
      keep <- integer()
    }
  }

  keep <- keep[order(tau[keep])]
  g <- g_modes[keep]
  if (plateau) g <- c(g, utils::tail(fit$g, 1L))

  list(
    g = g,
    tau = tau[keep],
    error = fit$error,
    keep = keep
  )
}

.experimentalPlacementWeights <- function(H, axis, s, wb, domain) {
  ns <- length(s)
  hs <- numeric(ns)
  hs[1L] <- 0.5 * log(s[2L] / s[1L])
  hs[ns] <- 0.5 * log(s[ns] / s[ns - 1L])
  hs[2:(ns - 1L)] <- 0.5 * (
    log(s[3:ns]) - log(s[1:(ns - 2L)])
  )

  kernel <- .experimentalMaxwellKernel(s, axis, domain)
  contribution <- sweep(kernel, 2, hs * exp(H), "*")
  prediction <- as.vector(kernel %*% (hs * exp(H)))

  rows <- seq_len(length(axis))
  contribution[rows, ] <- sweep(
    contribution[rows, , drop = FALSE],
    1, prediction[rows], "/"
  )

  wt <- colSums(contribution)
  wt <- wt / pracma::trapz(log(s), wt)
  (1 - wb) * wt + wb * mean(wt)
}

.experimentalCubic <- function(x, y, query) {
  n <- length(x)

  stopifnot(
    n >= 4L,
    length(y) == n,
    all(is.finite(x)),
    all(is.finite(y)),
    all(diff(x) > 0),
    all(is.finite(query)),
    all(query >= x[1L] & query <= x[n])
  )

  h <- diff(x)
  slopes <- diff(y) / h
  A <- matrix(0, n, n)
  rhs <- numeric(n)

  A[1L, 1:3] <- c(-h[2L], h[1L] + h[2L], -h[1L])
  A[n, (n - 2L):n] <- c(
    h[n - 1L],
    -(h[n - 2L] + h[n - 1L]),
    h[n - 2L]
  )

  for (i in 2:(n - 1L)) {
    A[i, (i - 1L):(i + 1L)] <- c(
      h[i - 1L], 2 * (h[i - 1L] + h[i]), h[i]
    )
    rhs[i] <- 6 * (slopes[i] - slopes[i - 1L])
  }

  second <- as.vector(solve(A, rhs))
  interval <- findInterval(query, x, all.inside = TRUE)
  width <- h[interval]
  a <- (x[interval + 1L] - query) / width
  b <- (query - x[interval]) / width

  a * y[interval] + b * y[interval + 1L] +
    ((a^3 - a) * second[interval] +
     (b^3 - b) * second[interval + 1L]) * width^2 / 6
}

.experimentalGridDensity <- function(x, px, N) {
  stopifnot(
    length(N) == 1L,
    N >= 1L,
    N == as.integer(N)
  )

  xi <- seq(min(x), max(x), length.out = 100L)
  pint <- .experimentalCubic(x, px, xi)
  ci <- as.vector(pracma::cumtrapz(xi, pint))

  stopifnot(is.finite(utils::tail(ci, 1L)), utils::tail(ci, 1L) > 0)
  ci <- ci / utils::tail(ci, 1L)
  stopifnot(all(diff(ci) > 0))

  if (N == 1L) {
    return(list(z = max(x), h = 0))
  }

  alpha <- 1 / (N - 1L)
  z <- numeric(N)
  z[1L] <- min(x)
  z[N] <- max(x)

  boundaries <- numeric(N + 1L)
  boundaries[1L] <- z[1L]
  boundaries[N + 1L] <- z[N]

  boundaries[2:N] <- .experimentalCubic(
    ci, xi, ((seq_len(N - 1L)) - 0.5) * alpha
  )

  if (N > 2L) {
    z[2:(N - 1L)] <- .experimentalCubic(
      ci, xi, seq_len(N - 2L) * alpha
    )
  }

  list(z = z, h = diff(boundaries))
}

.experimentalJacobian2Point <- function(fun, x, lb, ub) {
  lb <- rep_len(lb, length(x))
  ub <- rep_len(ub, length(x))
  f0 <- fun(x)

  h <- sqrt(.Machine$double.eps) *
    ifelse(x >= 0, 1, -1) * pmax(1, abs(x))

  lower_distance <- x - lb
  upper_distance <- ub - x
  violated <- x + h < lb | x + h > ub
  fitting <- abs(h) <= pmax(lower_distance, upper_distance)

  h[violated & fitting] <- -h[violated & fitting]

  forward <- !fitting & upper_distance >= lower_distance
  backward <- !fitting & upper_distance < lower_distance
  h[forward] <- upper_distance[forward]
  h[backward] <- -lower_distance[backward]

  J <- matrix(0, length(f0), length(x))

  for (j in seq_along(x)) {
    xp <- x
    xp[j] <- x[j] + h[j]
    dx <- xp[j] - x[j]

    if (!is.finite(dx) || dx == 0) {
      stop("Numerical Jacobian step is zero or non-finite.")
    }

    J[, j] <- (fun(xp) - f0) / dx
  }

  J
}

.experimentalPopulationSD <- function(x) {
  sqrt(mean((x - mean(x))^2))
}
