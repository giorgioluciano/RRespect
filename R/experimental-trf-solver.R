.trfCheckTermination <- function(
  dF, F, dx_norm, x_norm, ratio, ftol, xtol
) {
  f_ok <- dF < ftol * F && ratio > 0.25
  x_ok <- dx_norm < xtol * (xtol + x_norm)

  if (f_ok && x_ok) {
    4L
  } else if (f_ok) {
    2L
  } else if (x_ok) {
    3L
  } else {
    0L
  }
}

.trfNoBounds <- function(
  fun,
  jac,
  x0,
  ftol = 1e-8,
  xtol = 1e-8,
  gtol = 1e-8,
  max_nfev = NULL
) {
  stopifnot(
    is.function(fun),
    is.function(jac),
    is.numeric(x0),
    length(x0) >= 1L,
    is.null(dim(x0)),
    all(is.finite(x0)),
    length(ftol) == 1L,
    length(xtol) == 1L,
    length(gtol) == 1L,
    all(is.finite(c(ftol, xtol, gtol))),
    all(c(ftol, xtol, gtol) >= 0),
    any(c(ftol, xtol, gtol) > .Machine$double.eps)
  )

  x <- as.double(x0)
  n <- length(x)

  if (is.null(max_nfev)) {
    max_nfev <- 100L * n
  }

  stopifnot(
    length(max_nfev) == 1L,
    is.finite(max_nfev),
    max_nfev >= 1,
    max_nfev == floor(max_nfev)
  )

  evaluate_fun <- function(point, expected_length = NULL) {
    value <- fun(point)

    if (
      !is.numeric(value) ||
      !is.null(dim(value)) ||
      length(value) == 0L ||
      (!is.null(expected_length) &&
       length(value) != expected_length)
    ) {
      stop("TRF residual vector has an invalid shape.", call. = FALSE)
    }

    as.double(value)
  }

  f <- evaluate_fun(x)
  m <- length(f)

  if (!all(is.finite(f))) {
    stop("TRF initial residuals are not finite.", call. = FALSE)
  }

  evaluate_jac <- function(point) {
    value <- jac(point)

    if (
      !is.matrix(value) ||
      !is.numeric(value) ||
      !identical(dim(value), c(as.integer(m), as.integer(n))) ||
      !all(is.finite(value))
    ) {
      stop("TRF Jacobian has an invalid shape or values.", call. = FALSE)
    }

    value
  }

  J <- evaluate_jac(x)
  nfev <- 1L
  njev <- 1L
  cost <- 0.5 * sum(f^2)

  if (!is.finite(cost)) {
    stop("TRF initial cost is not finite.", call. = FALSE)
  }

  norm2 <- function(value) sqrt(sum(value^2))
  g <- as.vector(crossprod(J, f))
  Delta <- norm2(x)

  if (!is.finite(Delta)) {
    stop("TRF initial radius is not finite.", call. = FALSE)
  }

  if (Delta == 0) {
    Delta <- 1
  }

  alpha <- 0
  status <- 0L
  iteration <- 0L
  cost_history <- cost

  repeat {
    optimality <- max(abs(g))

    if (optimality < gtol) {
      status <- 1L
    }

    if (status != 0L || nfev >= max_nfev) {
      break
    }

    decomp <- svd(J)
    uf <- as.vector(crossprod(decomp$u, f))
    actual_reduction <- -1

    while (actual_reduction <= 0 && nfev < max_nfev) {
      if (!is.finite(Delta) || Delta <= 0) {
        stop("TRF radius collapsed numerically.", call. = FALSE)
      }

      subproblem <- .trfSolveLSQTrustRegion(
        n = n,
        m = m,
        uf = uf,
        s = decomp$d,
        V = decomp$v,
        Delta = Delta,
        initial_alpha = alpha
      )

      step <- subproblem$p
      alpha <- subproblem$alpha

      predicted_reduction <- -.trfEvaluateQuadratic(
        J, g, step
      )

      x_new <- x + step
      f_new <- evaluate_fun(x_new, m)
      nfev <- nfev + 1L
      step_norm <- norm2(step)

      if (!all(is.finite(f_new))) {
        Delta <- 0.25 * step_norm
        next
      }

      cost_new <- 0.5 * sum(f_new^2)

      if (!is.finite(cost_new)) {
        Delta <- 0.25 * step_norm
        next
      }

      actual_reduction <- cost - cost_new

      radius <- .trfUpdateRadius(
        Delta,
        actual_reduction,
        predicted_reduction,
        step_norm,
        step_norm > 0.95 * Delta
      )

      status <- .trfCheckTermination(
        actual_reduction,
        cost,
        step_norm,
        norm2(x),
        radius$ratio,
        ftol,
        xtol
      )

      if (status != 0L) {
        break
      }

      if (radius$Delta <= 0) {
        stop("TRF radius collapsed numerically.", call. = FALSE)
      }

      alpha <- alpha * Delta / radius$Delta
      Delta <- radius$Delta
    }

    if (actual_reduction > 0) {
      x <- x_new
      f <- f_new
      cost <- cost_new
      J <- evaluate_jac(x)
      njev <- njev + 1L
      g <- as.vector(crossprod(J, f))
      cost_history <- c(cost_history, cost)
    }

    iteration <- iteration + 1L
  }

  list(
    x = x,
    fun = f,
    jac = J,
    grad = g,
    cost = cost,
    optimality = max(abs(g)),
    nfev = nfev,
    njev = njev,
    nit = iteration,
    status = status,
    success = status > 0L,
    cost_history = cost_history
  )
}