.trfCLScaling <- function(x, g, lb, ub) {
  n <- length(x)

  stopifnot(
    length(g) == n,
    length(lb) == n,
    length(ub) == n,
    all(is.finite(x)),
    all(is.finite(g)),
    !anyNA(lb),
    !anyNA(ub),
    all(lb < ub),
    all(x >= lb & x <= ub)
  )

  v <- rep(1, n)
  dv <- rep(0, n)

  upper <- g < 0 & is.finite(ub)
  v[upper] <- ub[upper] - x[upper]
  dv[upper] <- -1

  lower <- g > 0 & is.finite(lb)
  v[lower] <- x[lower] - lb[lower]
  dv[lower] <- 1

  list(v = v, dv = dv)
}

.trfStepSizeToBound <- function(x, s, lb, ub) {
  n <- length(x)

  stopifnot(
    length(s) == n,
    length(lb) == n,
    length(ub) == n,
    all(is.finite(x)),
    all(is.finite(s)),
    !anyNA(lb),
    !anyNA(ub),
    all(lb < ub)
  )

  steps <- rep(Inf, n)
  moving <- s != 0

  steps[moving] <- pmax(
    (lb[moving] - x[moving]) / s[moving],
    (ub[moving] - x[moving]) / s[moving]
  )

  step <- min(steps)
  hits <- as.integer(steps == step) * as.integer(sign(s))

  list(step = step, hits = hits)
}

.trfUpdateRadius <- function(
    Delta,
    actual_reduction,
    predicted_reduction,
    step_norm,
    bound_hit
) {
  stopifnot(
    length(Delta) == 1L,
    is.finite(Delta),
    Delta > 0,
    length(actual_reduction) == 1L,
    is.finite(actual_reduction),
    length(predicted_reduction) == 1L,
    is.finite(predicted_reduction),
    length(step_norm) == 1L,
    is.finite(step_norm),
    step_norm >= 0,
    is.logical(bound_hit),
    length(bound_hit) == 1L,
    !is.na(bound_hit)
  )

  ratio <- if (predicted_reduction > 0) {
    actual_reduction / predicted_reduction
  } else if (
    predicted_reduction == 0 &&
    actual_reduction == 0
  ) {
    1
  } else {
    0
  }

  if (ratio < 0.25) {
    Delta <- 0.25 * step_norm
  } else if (ratio > 0.75 && bound_hit) {
    Delta <- 2 * Delta
  }

  list(Delta = Delta, ratio = ratio)
}


.trfSolveLSQTrustRegion <- function(
    n,
    m,
    uf,
    s,
    V,
    Delta,
    initial_alpha = NULL,
    rtol = 0.01,
    max_iter = 10L
) {
  stopifnot(
    n >= 1L,
    m >= 1L,
    length(s) %in% c(min(m, n), n),
    length(uf) == length(s),
    is.matrix(V),
    nrow(V) == n,
    ncol(V) == length(s),
    all(is.finite(s)),
    all(s >= 0),
    all(is.finite(uf)),
    all(is.finite(V)),
    length(Delta) == 1L,
    is.finite(Delta),
    Delta > 0,
    length(rtol) == 1L,
    is.finite(rtol),
    rtol > 0,
    length(max_iter) == 1L,
    is.finite(max_iter),
    max_iter >= 1L,
    max_iter == as.integer(max_iter)
  )

  if (!is.null(initial_alpha)) {
    stopifnot(
      length(initial_alpha) == 1L,
      is.finite(initial_alpha)
    )
  }

  norm2 <- function(x) sqrt(sum(x * x))

  suf <- s * uf
  full_rank <- m >= n &&
    s[length(s)] > .Machine$double.eps * m * s[1L]

  if (full_rank) {
    p <- -as.vector(V %*% (uf / s))

    if (norm2(p) <= Delta) {
      return(list(p = p, alpha = 0, n_iter = 0L))
    }
  }

  # A stationary quadratic model needs no step.
  if (all(suf == 0)) {
    return(list(
      p = rep(0, n),
      alpha = 0,
      n_iter = 0L
    ))
  }

  phi_and_derivative <- function(alpha) {
    denom <- s^2 + alpha
    p_norm <- norm2(suf / denom)

    list(
      phi = p_norm - Delta,
      derivative = -sum(suf^2 / denom^3) / p_norm
    )
  }

  alpha_upper <- norm2(suf) / Delta

  if (full_rank) {
    pd <- phi_and_derivative(0)
    alpha_lower <- -pd$phi / pd$derivative
  } else {
    alpha_lower <- 0
  }

  choose_alpha <- function(lower, upper) {
    max(0.001 * upper, sqrt(lower * upper))
  }

  if (
    is.null(initial_alpha) ||
    (!full_rank && initial_alpha == 0)
  ) {
    alpha <- choose_alpha(alpha_lower, alpha_upper)
  } else {
    alpha <- initial_alpha
  }

  for (it in seq_len(max_iter)) {
    if (alpha < alpha_lower || alpha > alpha_upper) {
      alpha <- choose_alpha(alpha_lower, alpha_upper)
    }

    pd <- phi_and_derivative(alpha)

    if (pd$phi < 0) {
      alpha_upper <- alpha
    }

    ratio <- pd$phi / pd$derivative
    alpha_lower <- max(alpha_lower, alpha - ratio)
    alpha <- alpha - (pd$phi + Delta) * ratio / Delta

    if (abs(pd$phi) < rtol * Delta) {
      break
    }
  }

  p <- -as.vector(V %*% (suf / (s^2 + alpha)))
  p_norm <- norm2(p)

  if (!all(is.finite(p)) || !is.finite(alpha) || p_norm == 0) {
    stop("TRF trust-region subproblem failed numerically.", call. = FALSE)
  }

  p <- p * (Delta / p_norm)

  list(
    p = p,
    alpha = alpha,
    n_iter = as.integer(it)
  )
}

.trfEvaluateQuadratic <- function(J, g, s, diag_h = NULL) {
  stopifnot(
    is.matrix(J),
    length(g) == ncol(J),
    length(s) == ncol(J),
    all(is.finite(J)),
    all(is.finite(g)),
    all(is.finite(s))
  )

  Js <- as.vector(J %*% s)
  q <- sum(Js^2)

  if (!is.null(diag_h)) {
    stopifnot(
      length(diag_h) == length(s),
      all(is.finite(diag_h))
    )
    q <- q + sum(diag_h * s^2)
  }

  0.5 * q + sum(g * s)
}

.trfBuildQuadratic1D <- function(
  J,
  g,
  s,
  diag_h = NULL,
  s0 = NULL
) {
  stopifnot(
    is.matrix(J),
    length(g) == ncol(J),
    length(s) == ncol(J),
    all(is.finite(J)),
    all(is.finite(g)),
    all(is.finite(s))
  )

  if (is.null(diag_h)) {
    diag_h <- rep(0, length(s))
  }

  stopifnot(
    length(diag_h) == length(s),
    all(is.finite(diag_h))
  )

  v <- as.vector(J %*% s)

  a <- 0.5 * (sum(v^2) + sum(diag_h * s^2))
  b <- sum(g * s)

  if (is.null(s0)) {
    return(list(a = a, b = b))
  }

  stopifnot(
    length(s0) == length(s),
    all(is.finite(s0))
  )

  u <- as.vector(J %*% s0)

  b <- b + sum(u * v) + sum(diag_h * s0 * s)
  c0 <- 0.5 * (sum(u^2) + sum(diag_h * s0^2)) +
    sum(g * s0)

  list(a = a, b = b, c = c0)
}

.trfMinimizeQuadratic1D <- function(a, b, lb, ub, c = 0) {
  stopifnot(
    length(a) == 1L,
    length(b) == 1L,
    length(c) == 1L,
    length(lb) == 1L,
    length(ub) == 1L,
    all(is.finite(c(a, b, c, lb, ub))),
    lb <= ub
  )

  candidates <- c(lb, ub)

  if (a != 0) {
    extremum <- -0.5 * b / a

    if (lb < extremum && extremum < ub) {
      candidates <- c(candidates, extremum)
    }
  }

  values <- candidates * (a * candidates + b) + c
  best <- which.min(values)

  list(t = candidates[best], value = values[best])
}

.trfIntersectTrustRegion <- function(x, s, Delta) {
  stopifnot(
    length(x) >= 1L,
    length(s) == length(x),
    all(is.finite(x)),
    all(is.finite(s)),
    length(Delta) == 1L,
    is.finite(Delta),
    Delta > 0
  )

  a <- sum(s^2)

  if (!is.finite(a) || a == 0) {
    stop(
      "Trust-region direction must have a finite nonzero norm.",
      call. = FALSE
    )
  }

  b <- sum(x * s)
  c0 <- sum(x^2) - Delta^2

  if (!is.finite(b) || !is.finite(c0)) {
    stop("Non-finite trust-region geometry.", call. = FALSE)
  }

  if (c0 > 0) {
    stop("Point is outside the trust region.", call. = FALSE)
  }

  discriminant <- b^2 - a * c0

  if (!is.finite(discriminant)) {
    stop("Non-finite trust-region discriminant.", call. = FALSE)
  }

  d <- sqrt(discriminant)
  signed_d <- if (b < 0) -d else d
  q <- -(b + signed_d)

  if (q == 0) {
    return(list(t_neg = 0, t_pos = 0))
  }

  roots <- sort(c(q / a, c0 / q))

  list(
    t_neg = roots[1L],
    t_pos = roots[2L]
  )
}

.trfSelectStep <- function(
  x, J_h, diag_h, g_h,
  p, p_h, d, Delta,
  lb, ub, theta
) {
  n <- length(x)

  stopifnot(
    n >= 1L,
    is.matrix(J_h),
    ncol(J_h) == n,
    length(diag_h) == n,
    length(g_h) == n,
    length(p) == n,
    length(p_h) == n,
    length(d) == n,
    length(lb) == n,
    length(ub) == n,
    all(is.finite(x)),
    all(is.finite(J_h)),
    all(is.finite(diag_h)),
    all(diag_h >= 0),
    all(is.finite(g_h)),
    all(is.finite(p)),
    all(is.finite(p_h)),
    all(is.finite(d)),
    all(d > 0),
    !anyNA(lb),
    !anyNA(ub),
    all(lb < ub),
    all(x >= lb & x <= ub),
    length(Delta) == 1L,
    is.finite(Delta),
    Delta > 0,
    length(theta) == 1L,
    is.finite(theta),
    theta > 0,
    theta < 1
  )

  quadratic <- function(step_h) {
    .trfEvaluateQuadratic(
      J_h, g_h, step_h, diag_h = diag_h
    )
  }

  result <- function(step, step_h, value, kind) {
    list(
      step = step,
      step_h = step_h,
      predicted_reduction = -value,
      kind = kind
    )
  }

  if (all(x + p >= lb & x + p <= ub)) {
    return(result(p, p_h, quadratic(p_h), "trust_region"))
  }

  bound <- .trfStepSizeToBound(x, p, lb, ub)
  p_stride <- bound$step

  r_h <- p_h
  hit <- bound$hits != 0L
  r_h[hit] <- -r_h[hit]
  r <- d * r_h

  p <- p * p_stride
  p_h <- p_h * p_stride
  x_on_bound <- x + p

  to_tr <- .trfIntersectTrustRegion(
    p_h, r_h, Delta
  )$t_pos

  to_bound <- .trfStepSizeToBound(
    x_on_bound, r, lb, ub
  )$step

  r_stride <- min(to_bound, to_tr)

  if (r_stride > 0) {
    r_lower <- (1 - theta) * p_stride / r_stride

    r_upper <- if (r_stride == to_bound) {
      theta * to_bound
    } else {
      to_tr
    }
  } else {
    r_lower <- 0
    r_upper <- -1
  }

  r_value <- Inf

  if (r_lower <= r_upper) {
    coef <- .trfBuildQuadratic1D(
      J_h, g_h, r_h,
      diag_h = diag_h,
      s0 = p_h
    )

    optimum <- .trfMinimizeQuadratic1D(
      coef$a, coef$b,
      r_lower, r_upper,
      c = coef$c
    )

    r_h <- p_h + optimum$t * r_h
    r <- d * r_h
    r_value <- optimum$value
  }

  p <- theta * p
  p_h <- theta * p_h
  p_value <- quadratic(p_h)

  ag_h <- -g_h
  ag <- d * ag_h
  ag_norm <- sqrt(sum(ag_h^2))

  if (ag_norm == 0) {
    ag <- rep(0, n)
    ag_h <- rep(0, n)
    ag_value <- 0
  } else {
    to_tr <- Delta / ag_norm
    to_bound <- .trfStepSizeToBound(
      x, ag, lb, ub
    )$step

    ag_upper <- if (to_bound < to_tr) {
      theta * to_bound
    } else {
      to_tr
    }

    coef <- .trfBuildQuadratic1D(
      J_h, g_h, ag_h,
      diag_h = diag_h
    )

    optimum <- .trfMinimizeQuadratic1D(
      coef$a, coef$b, 0, ag_upper
    )

    ag_h <- optimum$t * ag_h
    ag <- optimum$t * ag
    ag_value <- optimum$value
  }

  if (p_value < r_value && p_value < ag_value) {
    result(p, p_h, p_value, "constrained")
  } else if (r_value < p_value && r_value < ag_value) {
    result(r, r_h, r_value, "reflected")
  } else {
    result(ag, ag_h, ag_value, "cauchy")
  }
}


.trfFindActiveConstraints <- function(x, lb, ub, rtol = 1e-10) {
  n <- length(x)

  stopifnot(
    n >= 1L,
    length(lb) == n,
    length(ub) == n,
    all(is.finite(x)),
    !anyNA(lb),
    !anyNA(ub),
    all(lb < ub),
    length(rtol) == 1L,
    is.finite(rtol),
    rtol >= 0
  )

  active <- integer(n)

  if (rtol == 0) {
    active[x <= lb] <- -1L
    active[x >= ub] <- 1L
    return(active)
  }

  lower_dist <- x - lb
  upper_dist <- ub - x

  lower_threshold <- rtol * pmax(1, abs(lb))
  upper_threshold <- rtol * pmax(1, abs(ub))

  lower_active <- is.finite(lb) &
    lower_dist <= pmin(upper_dist, lower_threshold)

  active[lower_active] <- -1L

  upper_active <- is.finite(ub) &
    upper_dist <= pmin(lower_dist, upper_threshold)

  active[upper_active] <- 1L

  active
}

.trfNextAfterFinite <- function(x, toward) {
  stopifnot(
    length(x) == length(toward),
    all(is.finite(x)),
    !anyNA(toward)
  )

  vapply(seq_along(x), function(i) {
    value <- x[i]
    target <- toward[i]

    if (value == target) {
      return(as.double(target))
    }

    if (value == 0) {
      bytes <- integer(8)
      bytes[1L] <- 1L

      if (target < 0) {
        bytes[8L] <- 128L
      }
    } else {
      bytes <- as.integer(writeBin(
        as.double(value),
        raw(),
        size = 8L,
        endian = "little"
      ))

      increment <- (target > value) == (value > 0)

      for (j in seq_len(8L)) {
        if (increment) {
          if (bytes[j] < 255L) {
            bytes[j] <- bytes[j] + 1L
            break
          }
          bytes[j] <- 0L
        } else {
          if (bytes[j] > 0L) {
            bytes[j] <- bytes[j] - 1L
            break
          }
          bytes[j] <- 255L
        }
      }
    }

    readBin(
      as.raw(bytes),
      what = "double",
      n = 1L,
      size = 8L,
      endian = "little"
    )
  }, numeric(1))
}

.trfMakeStrictlyFeasible <- function(
  x, lb, ub, rstep = 1e-10
) {
  active <- .trfFindActiveConstraints(
    x, lb, ub, rtol = rstep
  )

  lower <- active == -1L
  upper <- active == 1L
  out <- as.double(x)

  if (rstep == 0) {
    out[lower] <- .trfNextAfterFinite(
      lb[lower], ub[lower]
    )
    out[upper] <- .trfNextAfterFinite(
      ub[upper], lb[upper]
    )
  } else {
    out[lower] <- lb[lower] +
      rstep * pmax(1, abs(lb[lower]))

    out[upper] <- ub[upper] -
      rstep * pmax(1, abs(ub[upper]))
  }

  tight <- out < lb | out > ub

  out[tight] <- 0.5 * (
    lb[tight] + ub[tight]
  )

  out
}
