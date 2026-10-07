
test_that("bound geometry tolerates an intermediate round-off violation", {
  lower <- 0.0002
  x <- lower - 2.25e-15

  result <- .trfStepSizeToBound(
    x = x,
    s = 1,
    lb = lower,
    ub = 500
  )

  expect_equal(result$step, 500 - x, tolerance = 1e-12)
  expect_identical(result$hits, 1L)
})

test_that("negative initial alpha is rebracketed for a deficient model", {
  J <- diag(c(2, 0))
  f <- c(2, 1)
  decomposition <- svd(J)

  result <- .trfSolveLSQTrustRegion(
    n = 2L,
    m = 2L,
    uf = as.vector(crossprod(decomposition$u, f)),
    s = decomposition$d,
    V = decomposition$v,
    Delta = 0.5,
    initial_alpha = -44348.774862900049
  )

  expect_equal(result$p, c(-0.5, 0), tolerance = 1e-10)
  expect_equal(result$alpha, 4, tolerance = 1e-8)
})

test_that("bounded TRF handles an inactive rank-deficient coordinate", {
  result <- .trfBounds(
    fun = function(x) c(x[1] - 1, 0),
    jac = function(x) diag(c(1, 0)),
    x0 = c(0.5, 0.5),
    lb = 0,
    ub = 2
  )

  expect_true(result$success)
  expect_equal(result$x, c(1, 0.5), tolerance = 1e-6)
  expect_lt(result$cost, 1e-14)
  expect_true(all(result$x > 0 & result$x < 2))
})
