
test_that("experimental error band uses selected snapshots", {
  H <- rbind(c(0, 2, 100), c(3, 3, 10))
  p <- c(0.45, 0.45, 0.1)

  expect_equal(.respectErrorBand(H, p), c(1, 0),
               tolerance = 1e-12)
})

test_that("experimental error band handles one selected snapshot", {
  H <- rbind(c(0, 2, 100), c(3, 3, 10))
  p <- c(0.9, 0.05, 0.05)

  expect_equal(.respectErrorBand(H, p), c(0, 0),
               tolerance = 1e-12)
})

test_that("experimental error band falls back to weighted variance", {
  x <- seq(0, 1, length.out = 50)
  H <- rbind(x, rep(0, 50))
  p <- rep(1 / 50, 50)

  expected <- c(
    sqrt(sum(p * (x - sum(p * x))^2)),
    0
  )

  out <- .respectErrorBand(H, p)
  expect_true(all(is.finite(out)))
  expect_equal(out, expected, tolerance = 1e-12)
})

test_that("lambda density defaults depend on the selected solver", {
  expect_identical(
    setParams(domain = "time", solver = "legacy")$lamDensity,
    3L
  )
  expect_identical(
    setParams(domain = "time", solver = "experimental")$lamDensity,
    2L
  )
  expect_identical(
    setParams(
      domain = "time",
      solver = "experimental",
      lamDensity = 5L
    )$lamDensity,
    5L
  )
})
