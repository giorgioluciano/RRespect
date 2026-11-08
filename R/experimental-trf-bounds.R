.trfBounds <-
function (fun, jac, x0, lb = -Inf, ub = Inf, ftol = 1e-08, xtol = 1e-08, gtol = 1e-08, max_nfev = NULL)
{
    stopifnot(is.function(fun), is.function(jac), is.numeric(x0), is.null(dim(x0)), length(x0) >= 1L,
        all(is.finite(x0)), length(ftol) == 1L, length(xtol) == 1L, length(gtol) == 1L, all(is.finite(c(ftol,
            xtol, gtol))), all(c(ftol, xtol, gtol) >= 0), any(c(ftol, xtol, gtol) > .Machine$double.eps))
    n <- length(x0)
    expand_bound <- function(value) {
        stopifnot(is.numeric(value), length(value) %in% c(1L, n), !anyNA(value))
        rep(as.double(value), length.out = n)
    }
    lb <- expand_bound(lb)
    ub <- expand_bound(ub)
    if (any(lb >= ub)) {
        stop("Each lower bound must be below its upper bound.", call. = FALSE)
    }
    if (any(x0 < lb | x0 > ub)) {
        stop("Initial point is outside the bounds.", call. = FALSE)
    }
    if (all(lb == -Inf) && all(ub == Inf)) {
        return(.trfNoBounds(fun, jac, x0, ftol = ftol, xtol = xtol, gtol = gtol, max_nfev = max_nfev))
    }
    x <- .trfMakeStrictlyFeasible(as.double(x0), lb, ub)
    if (any(x <= lb | x >= ub)) {
        stop("No strictly interior starting point was obtained.", call. = FALSE)
    }
    if (is.null(max_nfev)) {
        max_nfev <- 100L * n
    }
    stopifnot(length(max_nfev) == 1L, is.finite(max_nfev), max_nfev >= 1, max_nfev == floor(max_nfev))
    evaluate_fun <- function(point, expected = NULL) {
        value <- fun(point)
        if (!is.numeric(value) || !is.null(dim(value)) || length(value) == 0L || (!is.null(expected) &&
            length(value) != expected)) {
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
        if (!is.matrix(value) || !is.numeric(value) || !identical(dim(value), c(as.integer(m), as.integer(n))) ||
            !all(is.finite(value))) {
            stop("TRF Jacobian has an invalid shape or values.", call. = FALSE)
        }
        value
    }
    norm2 <- function(value) sqrt(sum(value^2))
    J <- evaluate_jac(x)
    g <- as.vector(crossprod(J, f))
    cost <- 0.5 * sum(f^2)
    if (!is.finite(cost) || !all(is.finite(g))) {
        stop("TRF initial cost or gradient is not finite.", call. = FALSE)
    }
    scaling <- .trfCLScaling(x, g, lb, ub)
    Delta <- norm2(x/sqrt(scaling$v))
    if (!is.finite(Delta)) {
        stop("TRF initial radius is not finite.", call. = FALSE)
    }
    if (Delta == 0) {
        Delta <- 1
    }
    nfev <- 1L
    njev <- 1L
    iteration <- 0L
    alpha <- 0
    status <- 0L
    cost_history <- cost
    repeat {
        scaling <- .trfCLScaling(x, g, lb, ub)
        optimality <- max(abs(g * scaling$v))
        if (optimality < gtol) {
            status <- 1L
        }
        if (status != 0L || nfev >= max_nfev) {
            break
        }
        d <- sqrt(scaling$v)
        diag_h <- g * scaling$dv
        g_h <- d * g
        J_h <- sweep(J, 2L, d, "*")
        J_augmented <- rbind(J_h, diag(sqrt(diag_h), nrow = n, ncol = n))
        f_augmented <- c(f, rep(0, n))
        decomp <- svd(J_augmented)
        uf <- as.vector(crossprod(decomp$u, f_augmented))
        theta <- max(0.995, 1 - optimality)
        actual_reduction <- -1
        while (actual_reduction <= 0 && nfev < max_nfev) {
            if (!is.finite(Delta) || Delta <= 0) {
                stop("TRF radius collapsed numerically.", call. = FALSE)
            }
            subproblem <- .trfSolveLSQTrustRegion(n = n, m = m, uf = uf, s = decomp$d, V = decomp$v,
                Delta = Delta, initial_alpha = alpha)
            alpha <- subproblem$alpha
            p_h <- subproblem$p
            selected <- .trfSelectStep(x, J_h, diag_h, g_h, p = d * p_h, p_h = p_h, d = d, Delta = Delta,
                lb = lb, ub = ub, theta = theta)
            x_new <- .trfMakeStrictlyFeasible(x + selected$step, lb, ub, rstep = 0)
            if (any(x_new <= lb | x_new >= ub)) {
                stop("TRF trial point is not strictly feasible.", call. = FALSE)
            }
            f_new <- evaluate_fun(x_new, m)
            nfev <- nfev + 1L
            step_h_norm <- norm2(selected$step_h)
            if (!all(is.finite(f_new))) {
                Delta <- 0.25 * step_h_norm
                next
            }
            cost_new <- 0.5 * sum(f_new^2)
            if (!is.finite(cost_new)) {
                Delta <- 0.25 * step_h_norm
                next
            }
            actual_reduction <- cost - cost_new
            radius <- .trfUpdateRadius(Delta, actual_reduction, selected$predicted_reduction, step_h_norm,
                step_h_norm > 0.95 * Delta)
            status <- .trfCheckTermination(actual_reduction, cost, norm2(selected$step), norm2(x), radius$ratio,
                ftol, xtol)
            if (status != 0L) {
                break
            }
            if (radius$Delta <= 0) {
                stop("TRF radius collapsed numerically.", call. = FALSE)
            }
            alpha <- alpha * Delta/radius$Delta
            Delta <- radius$Delta
        }
        if (actual_reduction > 0) {
            x <- x_new
            f <- f_new
            cost <- cost_new
            J <- evaluate_jac(x)
            njev <- njev + 1L
            g <- as.vector(crossprod(J, f))
            if (!all(is.finite(g))) {
                stop("TRF gradient is not finite.", call. = FALSE)
            }
            cost_history <- c(cost_history, cost)
        }
        iteration <- iteration + 1L
    }
    scaling <- .trfCLScaling(x, g, lb, ub)
    list(x = x, fun = f, jac = J, grad = g, cost = cost, optimality = max(abs(g * scaling$v)), active_mask = .trfFindActiveConstraints(x,
        lb, ub, rtol = xtol), nfev = nfev, njev = njev, nit = iteration, status = status, success = status >
        0L, cost_history = cost_history)
}
