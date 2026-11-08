
test_that("experimental cubic reproduces a cubic polynomial", {
  x <- c(-2, -0.7, 0, 0.4, 1.8, 3)
  query <- seq(min(x), max(x), length.out = 41)
  f <- function(x) x^3 - 2 * x + 1

  expect_equal(
    .experimentalCubic(x, f(x), query),
    f(query),
    tolerance = 1e-11
  )
})

test_that("experimental density grid handles one and two nodes", {
  x <- seq(-2, 3, length.out = 10)
  px <- rep(1, length(x))

  one <- .experimentalGridDensity(x, px, 1L)
  two <- .experimentalGridDensity(x, px, 2L)
  five <- .experimentalGridDensity(x, px, 5L)

  expect_equal(one$z, 3)
  expect_equal(one$h, 0)
  expect_equal(two$z, c(-2, 3))
  expect_equal(two$h, c(2.5, 2.5), tolerance = 1e-10)
  expect_equal(five$z, seq(-2, 3, length.out = 5),
               tolerance = 1e-10)
})

test_that("experimental population SD uses the population denominator", {
  x <- c(1, 2, 3, 4)
  expect_equal(.experimentalPopulationSD(x), sqrt(1.25))
  expect_equal(
    .experimentalPopulationSD(x),
    sd(x) * sqrt(3 / 4)
  )
})

test_that("experimental numerical Jacobian respects active bounds", {
  A <- rbind(c(1, 2), c(3, -1), c(0.5, 4))
  fun <- function(x) as.vector(A %*% x)

  for (x in list(c(0, 1), c(0.4, 0.6))) {
    J <- .experimentalJacobian2Point(fun, x, lb = 0, ub = 1)
    expect_equal(J, A, tolerance = 1e-7)
  }
})

test_that("experimental NNLS recovers known modes and plateau", {
  axis <- exp(seq(log(0.01), log(100), length.out = 40))
  tau <- c(0.5, 10)
  g <- c(1, 0.4)

  for (domain in c("time", "frequency")) {
    K <- .experimentalMaxwellKernel(tau, axis, domain)

    for (plateau in c(FALSE, TRUE)) {
      G <- as.vector(K %*% g)

      if (plateau) {
        if (domain == "time") {
          G <- G + 0.2
        } else {
          G[seq_along(axis)] <- G[seq_along(axis)] + 0.2
        }
      }

      fit <- .experimentalNnLS(
        tau, axis, G, rep(1, length(G)), plateau, domain
      )

      expect_equal(
        fit$g,
        if (plateau) c(g, 0.2) else g,
        tolerance = 1e-8
      )
      expect_lt(fit$error, 1e-14)
    }
  }
})

test_that("experimental tau windows depend on the domain", {
  expect_equal(
    .experimentalTauBounds(c(0.1, 100), "time"),
    c(0.002, 5000)
  )
  expect_equal(
    .experimentalTauBounds(c(0.1, 100), "frequency"),
    c(0.0002, 500)
  )
})
