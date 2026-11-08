
test_that("bounded TRF finds an interior scalar solution", {
  out <- .trfBounds(
    fun = function(x) x - 0.3,
    jac = function(x) matrix(1),
    x0 = 0.8, lb = 0, ub = 1
  )

  expect_true(out$success)
  expect_equal(out$x, 0.3, tolerance = 1e-6)
  expect_lt(out$cost, 1e-12)
  expect_true(all(diff(out$cost_history) <= 0))
})

test_that("bounded TRF approaches an upper-bound optimum", {
  out <- .trfBounds(
    fun = function(x) x - 2,
    jac = function(x) matrix(1),
    x0 = 0.5, lb = 0, ub = 1
  )

  expect_true(out$success)
  expect_equal(out$x, 1, tolerance = 1e-6)
  expect_true(out$x > 0 && out$x < 1)
  expect_lt(out$optimality, 1e-7)
})

test_that("bounded TRF handles mixed bounds and scalar expansion", {
  target <- c(-1, 2)

  out <- .trfBounds(
    fun = function(x) x - target,
    jac = function(x) diag(2),
    x0 = c(0.5, 0.5),
    lb = 0, ub = c(Inf, 1)
  )

  expect_true(out$success)
  expect_equal(out$x, c(0, 1), tolerance = 1e-6)
  expect_true(all(out$x > 0))
  expect_lt(out$x[2], 1)
})

test_that("bounded TRF moves a bound start inward with ftol disabled", {
  out <- .trfBounds(
    fun = function(x) x - 0.5,
    jac = function(x) matrix(1),
    x0 = 0, lb = 0, ub = 1,
    ftol = 0,
    max_nfev = 1000L
  )

  expect_true(out$success)
  expect_equal(out$x, 0.5, tolerance = 1e-6)
})

test_that("bounded TRF rejects infeasible starts and invalid intervals", {
  expect_error(
    .trfBounds(
      function(x) x,
      function(x) matrix(1),
      x0 = 2, lb = 0, ub = 1
    ),
    "outside the bounds"
  )

  expect_error(
    .trfBounds(
      function(x) x,
      function(x) matrix(1),
      x0 = 1, lb = 1, ub = 1
    ),
    "lower bound"
  )
})

test_that("bounded TRF evaluation limit is not a success", {
  out <- .trfBounds(
    fun = function(x) x - 2,
    jac = function(x) matrix(1),
    x0 = 0.5, lb = 0, ub = 1,
    max_nfev = 1L
  )

  expect_false(out$success)
  expect_identical(out$status, 0L)
  expect_identical(out$nfev, 1L)
})
