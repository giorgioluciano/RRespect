test_that("Coleman-Li scaling follows gradient and finite bounds", {
    out <- .trfCLScaling(x = c(2, 2, 2, 2), g = c(-3, 3, 0, -3), lb = c(0, 0, 0, -Inf), ub = c(5, 5,
        5, Inf))
    expect_equal(out$v, c(3, 2, 1, 1))
    expect_equal(out$dv, c(-1, 1, 0, 0))
})

test_that("Coleman-Li scaling handles active bounds", {
    out <- .trfCLScaling(x = c(0, 5), g = c(2, -2), lb = c(0, 0), ub = c(5, 5))
    expect_equal(out$v, c(0, 0))
    expect_equal(out$dv, c(1, -1))
})

test_that("first bound hit has the correct stride and direction", {
    out <- .trfStepSizeToBound(x = c(1, 2), s = c(2, -1), lb = c(0, 0), ub = c(3, 5))
    expect_equal(out$step, 1)
    expect_identical(out$hits, c(1L, 0L))
})

test_that("simultaneous bound hits are retained", {
    out <- .trfStepSizeToBound(x = c(1, 2), s = c(2, -2), lb = c(0, 0), ub = c(3, 5))
    expect_equal(out$step, 1)
    expect_identical(out$hits, c(1L, -1L))
})

test_that("unbounded and stationary directions have infinite stride", {
    for (s in list(c(1, -1), c(0, 0))) {
        out <- .trfStepSizeToBound(x = c(1, 2), s = s, lb = c(-Inf, -Inf), ub = c(Inf, Inf))
        expect_identical(out$step, Inf)
    }
    stationary <- .trfStepSizeToBound(x = c(1, 2), s = c(0, 0), lb = c(0, 0), ub = c(3, 5))
    expect_identical(stationary$hits, c(0L, 0L))
})

test_that("trust-region radius shrinks after a poor prediction", {
    out <- .trfUpdateRadius(2, 0.1, 1, 0.8, FALSE)
    expect_equal(out$ratio, 0.1)
    expect_equal(out$Delta, 0.2)
})

test_that("trust-region radius grows only on a boundary-reaching step", {
    grow <- .trfUpdateRadius(2, 0.9, 1, 2, TRUE)
    keep <- .trfUpdateRadius(2, 0.9, 1, 1, FALSE)
    expect_equal(grow$Delta, 4)
    expect_equal(keep$Delta, 2)
})

test_that("zero reductions produce unit ratio", {
    out <- .trfUpdateRadius(2, 0, 0, 1, FALSE)
    expect_equal(out$ratio, 1)
    expect_equal(out$Delta, 2)
})

test_that("SVD subproblem accepts an interior Gauss-Newton step", {
    J <- diag(c(2, 1))
    f <- c(2, 2)
    decomp <- svd(J)
    out <- .trfSolveLSQTrustRegion(n = ncol(J), m = nrow(J), uf = as.vector(crossprod(decomp$u, f)),
        s = decomp$d, V = decomp$v, Delta = 3)
    expect_equal(out$p, c(-1, -2), tolerance = 1e-12)
    expect_equal(out$alpha, 0)
    expect_identical(out$n_iter, 0L)
})

test_that("SVD subproblem solves an isotropic boundary step", {
    J <- diag(2)
    f <- c(3, 4)
    decomp <- svd(J)
    out <- .trfSolveLSQTrustRegion(n = ncol(J), m = nrow(J), uf = as.vector(crossprod(decomp$u, f)),
        s = decomp$d, V = decomp$v, Delta = 1)
    expect_equal(out$p, c(-0.6, -0.8), tolerance = 1e-10)
    expect_equal(out$alpha, 4, tolerance = 1e-08)
    expect_equal(sqrt(sum(out$p^2)), 1, tolerance = 1e-12)
    expect_gt(out$n_iter, 0L)
})

test_that("SVD subproblem handles a rank-deficient Jacobian", {
    J <- diag(c(2, 0))
    f <- c(2, 1)
    decomp <- svd(J)
    out <- .trfSolveLSQTrustRegion(n = ncol(J), m = nrow(J), uf = as.vector(crossprod(decomp$u, f)),
        s = decomp$d, V = decomp$v, Delta = 0.5)
    expect_true(all(is.finite(out$p)))
    expect_equal(out$p, c(-0.5, 0), tolerance = 1e-10)
    expect_equal(out$alpha, 4, tolerance = 1e-08)
})

test_that("SVD subproblem returns zero for a stationary model", {
    out <- .trfSolveLSQTrustRegion(n = 2L, m = 2L, uf = c(0, 1), s = c(2, 0), V = diag(2), Delta = 1)
    expect_equal(out$p, c(0, 0))
    expect_equal(out$alpha, 0)
    expect_identical(out$n_iter, 0L)
})

test_that("quadratic evaluation includes the diagonal correction", {
    J <- diag(c(2, 1))
    g <- c(-1, 2)
    s <- c(1, -2)
    expect_equal(.trfEvaluateQuadratic(J, g, s), -1)
    expect_equal(.trfEvaluateQuadratic(J, g, s, diag_h = c(2, 3)), 6)
})

test_that("one-dimensional coefficients match direct evaluation", {
    J <- matrix(c(1, 2, -1, 3, 0, 2), nrow = 3)
    g <- c(-2, 1)
    s <- c(0.5, -1)
    s0 <- c(1, 0.25)
    diag_h <- c(0.2, 0.4)
    coef <- .trfBuildQuadratic1D(J, g, s, diag_h = diag_h, s0 = s0)
    for (t in c(-1, 0, 0.3, 2)) {
        polynomial <- t * (coef$a * t + coef$b) + coef$c
        direct <- .trfEvaluateQuadratic(J, g, s0 + t * s, diag_h = diag_h)
        expect_equal(polynomial, direct, tolerance = 1e-12)
    }
})

test_that("quadratic line through zero has no constant term", {
    coef <- .trfBuildQuadratic1D(J = diag(2), g = c(-3, 1), s = c(1, 2))
    expect_equal(coef$a, 2.5)
    expect_equal(coef$b, -1)
    expect_null(coef$c)
})

test_that("quadratic minimization chooses an interior minimum", {
    out <- .trfMinimizeQuadratic1D(a = 2, b = -4, lb = 0, ub = 3, c = 5)
    expect_equal(out$t, 1)
    expect_equal(out$value, 3)
})

test_that("quadratic minimization respects interval endpoints", {
    convex <- .trfMinimizeQuadratic1D(1, -4, 0, 1)
    linear <- .trfMinimizeQuadratic1D(0, -2, 0, 3)
    concave <- .trfMinimizeQuadratic1D(-1, 0, -1, 2)
    expect_equal(convex$t, 1)
    expect_equal(linear$t, 3)
    expect_equal(concave$t, 2)
})

test_that("quadratic minimization accepts a zero-width interval", {
    out <- .trfMinimizeQuadratic1D(1, -2, 0.5, 0.5)
    expect_equal(out$t, 0.5)
    expect_equal(out$value, -0.75)
})

test_that("trust-region intersection works from the origin", {
    out <- .trfIntersectTrustRegion(x = c(0, 0), s = c(3, 4), Delta = 2)
    expect_equal(out$t_neg, -0.4)
    expect_equal(out$t_pos, 0.4)
})

test_that("trust-region roots reach the sphere from an interior point", {
    x <- c(0.2, -0.3)
    s <- c(1, 2)
    Delta <- 1
    out <- .trfIntersectTrustRegion(x, s, Delta)
    expect_lt(out$t_neg, 0)
    expect_gt(out$t_pos, 0)
    for (t in c(out$t_neg, out$t_pos)) {
        expect_equal(sqrt(sum((x + t * s)^2)), Delta, tolerance = 1e-12)
    }
})

test_that("boundary intersections distinguish inward and outward directions", {
    inward <- .trfIntersectTrustRegion(x = c(1, 0), s = c(-1, 0), Delta = 1)
    outward <- .trfIntersectTrustRegion(x = c(1, 0), s = c(1, 0), Delta = 1)
    expect_equal(inward$t_neg, 0)
    expect_equal(inward$t_pos, 2)
    expect_equal(outward$t_neg, -2)
    expect_equal(outward$t_pos, 0)
})

test_that("a tangent boundary direction has a double zero root", {
    out <- .trfIntersectTrustRegion(x = c(1, 0), s = c(0, 1), Delta = 1)
    expect_equal(out$t_neg, 0)
    expect_equal(out$t_pos, 0)
})

test_that("invalid trust-region geometry is rejected", {
    expect_error(.trfIntersectTrustRegion(c(0, 0), c(0, 0), 1), "nonzero norm")
    expect_error(.trfIntersectTrustRegion(c(2, 0), c(1, 0), 1), "outside the trust region")
})

test_that("step selection preserves a feasible trust-region step", {
    out <- .trfSelectStep(x = c(1, 1), J_h = diag(2), diag_h = c(0, 0), g_h = c(-1, -1), p = c(0.2, 0.3),
        p_h = c(0.2, 0.3), d = c(1, 1), Delta = 1, lb = c(0, 0), ub = c(2, 2), theta = 0.995)
    expect_identical(out$kind, "trust_region")
    expect_equal(out$step, c(0.2, 0.3))
    expect_equal(out$step_h, c(0.2, 0.3))
    expect_equal(out$predicted_reduction, 0.435)
})

test_that("step selection matches the scalar Cauchy reference", {
    out <- .trfSelectStep(x = 0.9, J_h = matrix(1), diag_h = 0, g_h = 1, p = 0.5, p_h = 0.5, d = 1, Delta = 1,
        lb = 0, ub = 1, theta = 0.995)
    expect_identical(out$kind, "cauchy")
    expect_equal(out$step, -0.8955, tolerance = 1e-12)
    expect_equal(out$step_h, out$step)
    expect_equal(out$predicted_reduction, 0.494539875, tolerance = 1e-12)
})

test_that("step selection matches the two-dimensional reflected reference", {
    out <- .trfSelectStep(x = c(0.9, 0), J_h = diag(2), diag_h = c(0, 0), g_h = c(-0.5, -0.5), p = c(0.5,
        0.5), p_h = c(0.5, 0.5), d = c(1, 1), Delta = 1, lb = c(0, -Inf), ub = c(1, Inf), theta = 0.995)
    expect_identical(out$kind, "reflected")
    expect_equal(out$step, c(0.0996428571428571, 0.100357142857143), tolerance = 1e-12)
    expect_equal(out$predicted_reduction, 0.0899998724489796, tolerance = 1e-12)
})

test_that("step selection prefers a validated constrained candidate", {
    out <- .trfSelectStep(x = c(0.9, 0.2), J_h = diag(c(1.6, 2.6)), diag_h = c(0, 0), g_h = c(-1, 0.7),
        p = c(0.4, -0.4), p_h = c(0.4, -0.4), d = c(1, 1), Delta = 1, lb = c(0, 0), ub = c(1, 1), theta = 0.995)
    expect_identical(out$kind, "constrained")
    expect_equal(out$step, c(0.0995, -0.0995), tolerance = 1e-12)
    expect_equal(out$predicted_reduction, 0.123014835, tolerance = 1e-12)
})

test_that("scaled step selection preserves geometry and model consistency", {
    x <- c(0.9, 0.2)
    d <- c(2, 0.5)
    p_h <- c(0.3, -0.4)
    J_h <- diag(c(1, 2))
    diag_h <- c(0.1, 0.2)
    g_h <- c(-0.2, 1)
    out <- .trfSelectStep(x = x, J_h = J_h, diag_h = diag_h, g_h = g_h, p = d * p_h, p_h = p_h, d = d,
        Delta = 1, lb = c(0, 0), ub = c(1, 1), theta = 0.995)
    expect_true(all(x + out$step >= 0))
    expect_true(all(x + out$step <= 1))
    expect_lte(sqrt(sum(out$step_h^2)), 1 + 1e-12)
    expect_equal(out$step, d * out$step_h)
    expect_equal(out$predicted_reduction, -.trfEvaluateQuadratic(J_h, g_h, out$step_h, diag_h = diag_h),
        tolerance = 1e-12)
})

test_that("active constraints distinguish lower upper and interior points", {
    out <- .trfFindActiveConstraints(x = c(0, 1, 0.5, 2), lb = c(0, 0, 0, -Inf), ub = c(1, 1, 1, Inf))
    expect_identical(out, c(-1L, 1L, 0L, 0L))
})

test_that("active constraints use relative bound thresholds", {
    out <- .trfFindActiveConstraints(x = c(100 + 5e-07, 100 + 2e-06, 200 - 1e-06), lb = c(100, 100, 100),
        ub = c(200, 200, 200), rtol = 1e-08)
    expect_identical(out, c(-1L, 0L, 1L))
})

test_that("zero tolerance detects exact bounds and violations only", {
    out <- .trfFindActiveConstraints(x = c(0, 1, 1e-12, -0.1, 1.1), lb = rep(0, 5), ub = rep(1, 5), rtol = 0)
    expect_identical(out, c(-1L, 1L, 0L, -1L, 1L))
})

test_that("closest bound wins when activity thresholds overlap", {
    out <- .trfFindActiveConstraints(x = c(0.2, 0.8, 0.5), lb = c(0, 0, 0), ub = c(1, 1, 1), rtol = 1)
    expect_identical(out, c(-1L, 1L, 1L))
})

test_that("infinite bounds are never classified as active", {
    out <- .trfFindActiveConstraints(x = c(0, 1, -1), lb = c(-Inf, 0, -Inf), ub = c(Inf, Inf, 0))
    expect_identical(out, c(0L, 0L, 0L))
})

test_that("active constraints reject invalid bounds and tolerance", {
    expect_error(.trfFindActiveConstraints(0, 1, 1))
    expect_error(.trfFindActiveConstraints(0, -1, 1, rtol = -1))
})

test_that("nextafter matches binary64 neighbors around one", {
    eps <- .Machine$double.eps
    out <- .trfNextAfterFinite(c(1, 1, -1, -1), c(Inf, -Inf, Inf, -Inf))
    expect_identical(out, c(1 + eps, 1 - eps/2, -1 + eps/2, -1 - eps))
})

test_that("nextafter handles zero and the smallest subnormal", {
    tiny <- .Machine$double.xmin * .Machine$double.eps
    expect_gt(tiny, 0)
    expect_identical(.trfNextAfterFinite(c(0, 0), c(1, -1)), c(tiny, -tiny))
    expect_equal(.trfNextAfterFinite(c(tiny, -tiny), c(0, 0)), c(0, 0), tolerance = 0)
})

test_that("nextafter preserves equal finite targets", {
    expect_identical(.trfNextAfterFinite(c(1, -2, 0), c(1, -2, 0)), c(1, -2, 0))
})

test_that("relative feasibility moves active bounds inward", {
    out <- .trfMakeStrictlyFeasible(x = c(0, 1, 0.5), lb = c(0, 0, 0), ub = c(1, 1, 1), rstep = 1e-04)
    expect_equal(out, c(1e-04, 1 - 1e-04, 0.5))
    expect_true(all(out > 0 & out < 1))
})

test_that("zero-step feasibility uses immediate representable neighbors", {
    tiny <- .Machine$double.xmin * .Machine$double.eps
    out <- .trfMakeStrictlyFeasible(x = c(0, 1), lb = c(0, 0), ub = c(1, 1), rstep = 0)
    expect_identical(out, c(tiny, 1 - .Machine$double.eps/2))
    expect_true(all(out > 0 & out < 1))
})

test_that("tight intervals use the midpoint after an excessive shift", {
    out <- .trfMakeStrictlyFeasible(x = 0, lb = 0, ub = 1e-12, rstep = 1e-10)
    expect_equal(out, 5e-13, tolerance = 0)
    expect_true(out > 0 && out < 1e-12)
})

test_that("feasibility leaves unbounded interior points unchanged", {
    expect_identical(.trfMakeStrictlyFeasible(x = c(-2, 3), lb = c(-Inf, -Inf), ub = c(Inf, Inf)), c(-2,
        3))
})
