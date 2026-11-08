
test_that("unbounded TRF solves a linear least-squares problem", {
  A <- rbind(
    c(1, 0),
    c(0, 2),
    c(1, 1)
  )
  target <- c(2, -1)
  b <- as.vector(A %*% target)

  out <- .trfNoBounds(
    fun = function(x) as.vector(A %*% x - b),
    jac = function(x) A,
    x0 = c(0, 0)
  )

  expect_true(out$success)
  expect_equal(out$x, target, tolerance = 1e-8)
  expect_lt(out$cost, 1e-16)
  expect_true(all(diff(out$cost_history) <= 0))
})

test_that("unbounded TRF solves the Rosenbrock residual problem", {
  fun <- function(x) {
    c(10 * (x[2] - x[1]^2), 1 - x[1])
  }

  jac <- function(x) {
    rbind(
      c(-20 * x[1], 10),
      c(-1, 0)
    )
  }

  out <- .trfNoBounds(
    fun, jac,
    x0 = c(-1.2, 1),
    max_nfev = 1000L
  )

  expect_true(out$success)
  expect_equal(out$x, c(1, 1), tolerance = 1e-6)
  expect_lt(out$cost, 1e-12)
  expect_true(all(diff(out$cost_history) <= 0))
})

test_that("evaluation limit is not reported as convergence", {
  out <- .trfNoBounds(
    fun = function(x) x - 2,
    jac = function(x) matrix(1),
    x0 = 0,
    max_nfev = 1L
  )

  expect_false(out$success)
  expect_identical(out$status, 0L)
  expect_identical(out$nfev, 1L)
  expect_equal(out$x, 0)
})

test_that("stationary initial point terminates without another evaluation", {
  out <- .trfNoBounds(
    fun = function(x) x - 2,
    jac = function(x) matrix(1),
    x0 = 2
  )

  expect_true(out$success)
  expect_identical(out$status, 1L)
  expect_identical(out$nfev, 1L)
})

test_that("unbounded TRF rejects invalid initial residuals and Jacobians", {
  expect_error(
    .trfNoBounds(
      fun = function(x) Inf,
      jac = function(x) matrix(1),
      x0 = 0
    ),
    "initial residuals"
  )

  expect_error(
    .trfNoBounds(
      fun = function(x) x - 1,
      jac = function(x) matrix(1, 2, 1),
      x0 = 0
    ),
    "Jacobian"
  )
})
